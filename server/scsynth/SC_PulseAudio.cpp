/*
    SuperCollider real time audio synthesis system
    Copyright (c) 2002 James McCartney. All rights reserved.
    http://www.audiosynth.com

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#include "SC_CoreAudio.h"
#include <stdarg.h>
#include "SC_Prototypes.h"
#include "SC_HiddenWorld.h"
#include "SC_WorldOptions.h"
#include "SC_Time.hpp"
#include <math.h>
#include <stdlib.h>

#include <pulse/pulseaudio.h>
#define SC_PA_USE_DLL

int32 server_timeseed() { return timeSeed(); }

#ifdef SC_PA_USE_DLL
#    include "SC_TimeDLL.hpp"
// =====================================================================
// Timing

static inline int64 sc_PAOSCTime() { return OSCTime(getTime()); }

static inline double sc_PAOSCTimeSeconds() { return (uint64)sc_PAOSCTime() * kOSCtoSecs; }

int64 oscTimeNow() { return sc_PAOSCTime(); }

void initializeScheduler() {}

#else // SC_PA_USE_DLL

int64 gOSCoffset = 0;

static inline int64 GetCurrentOSCTime() { return OSCTime(getTime()); }

int64 oscTimeNow() { return GetCurrentOSCTime(); }

int64 PaStreamTimeToOSC(PaTime pa_time) {
    uint64 s, f;
    s = (uint64)pa_time;
    f = (uint64)((pa_time - s) * 1000000 * kMicrosToOSCunits);

    return (s << 32) + f;
}

void initializeScheduler() { gOSCoffset = GetCurrentOSCTime(); }

#endif // SC_PA_USE_DLL

struct sinkInfo
{
    std::string name;
    std::string description;
    std::string driver;
    int nChannels;
};




class SC_PulseAudioDriver : public SC_AudioDriver {
    int mInputChannelCount, mOutputChannelCount;
    // LLP: In pulseaudio, these are setup with different functions
    pa_stream* mStreamIn;
    pa_stream* mStreamOut;
    int64 mPaStreamStartupTime;
    int64 mPaStreamStartupTimeOSC;
#ifdef SC_PA_USE_DLL
    double mMaxOutputLatency;
    SC_TimeDLL mDLL;
#endif
    static void pa_state_cb(pa_context *c, void *userdata);
    static void sinkListCB(pa_context *c, const pa_sink_info *i, int eol, void *userdata);
    static void PulseAudioCallbackStatic(pa_stream *s, size_t length, void *userdata);
    void PulseAudioCallback(pa_stream *s, size_t length);
    pa_threaded_mainloop *pa_ml;
    pa_mainloop_api *pa_mlapi;
    pa_context *pa_ctx;
    pa_sample_spec ss;
    bool sinkListReady;
    int pa_ready;
    std::map<std::string, sinkInfo> sinkInfos;
    std::string outSinkToUse;
protected:
    // Driver interface methods
    virtual bool DriverSetup(int* outNumSamplesPerCallback, double* outSampleRate) override;
    virtual bool DriverStart() override;
    virtual bool DriverStop() override;

public:
    SC_PulseAudioDriver(struct World* inWorld);
    virtual ~SC_PulseAudioDriver();
};

SC_AudioDriver* SC_NewAudioDriver(struct World* inWorld) { 
    printf("In SC_NewAudioDriver\n");
    return new SC_PulseAudioDriver(inWorld); 
}

#define PRINT_PULSEAUDIO_ERROR(function, errorcode)                                                                     \
    printf("SC_PulseAudioDriver: PulseAudio failed at %s with error: '%s'\n", #function, pa_strerror(errorcode))

SC_PulseAudioDriver::SC_PulseAudioDriver(struct World* inWorld):
    SC_AudioDriver(inWorld),
    mStreamIn(0), mStreamOut(0), sinkListReady(false), pa_ready(0)
#ifdef SC_PA_USE_DLL
    ,
    mMaxOutputLatency(0.)
#endif
{
    printf("In SC_PulseAudioDriver::SC_PulseAudioDriver\n");
}

SC_PulseAudioDriver::~SC_PulseAudioDriver() {
    // LLP: the streams are closed at DriverStop time, so just need to terminate PulseAudio
    printf("In SC_PulseAudioDriver::~SC_PulseAudioDriver\n");
    pa_threaded_mainloop_stop(pa_ml);
    pa_threaded_mainloop_free(pa_ml);
}


#if 0
static void SC_PulseAudioStreamInCallback(pa_stream *p, size_t nbytes, void *userdata)
{
    SC_PulseAudioDriver* driver = (SC_PulseAudioDriver*)userData;

    return driver->PulseAudioCallback(input, output, frameCount, timeInfo, statusFlags);
}
#endif
void sc_SetDenormalFlags();


void SC_PulseAudioDriver::PulseAudioCallbackStatic(pa_stream *s, size_t length, void *userdata) 
{
    SC_PulseAudioDriver* driver = (SC_PulseAudioDriver*)userdata;

    driver->PulseAudioCallback(s, length);
}
void SC_PulseAudioDriver::PulseAudioCallback(pa_stream *s, size_t length)
{
    //printf("PulseAudioCallback!!\n");
    sc_SetDenormalFlags();
    World* world = mWorld;
#ifdef SC_PA_USE_DLL
    mDLL.Update(sc_PAOSCTimeSeconds());

#    if SC_PA_DEBUG_DLL
    static int tick = 0;
    if (++tick >= 10) {
        tick = 0;
        scprintf("DLL: t %.6f p %.9f sr %.6f e %.9f avg(e) %.9f inc %.9f\n", mDLL.PeriodTime(), mDLL.Period(),
                 mDLL.SampleRate(), mDLL.Error(), mDLL.AvgError(), mOSCincrement * kOSCtoSecs);

        // scprintf("mOSCbuftime1 %llu \t %llu \t %f \n",GetCurrentOSCTime(),(uint64)((mDLL.PeriodTime() -
        // mMaxOutputLatency) * kSecondsToOSCunits + .5),((mDLL.PeriodTime() - mMaxOutputLatency) * kSecondsToOSCunits +
        // .5));
    }
#    endif
#endif

    try {
#if !defined(SC_PA_USE_DLL)
        // synchronise against the output buffer - timeInfo->currentTime is 0.0 bug in PA?
        if (mPaStreamStartupTime == 0 && mPaStreamStartupTimeOSC == 0) {
            mPaStreamStartupTimeOSC = GetCurrentOSCTime();
            mPaStreamStartupTime = timeInfo->outputBufferDacTime;
        }
        mOSCbuftime = PaStreamTimeToOSC(timeInfo->outputBufferDacTime - mPaStreamStartupTime) + mPaStreamStartupTimeOSC;
#endif

        mFromEngine.Free();
        mToEngine.Perform();
        mOscPacketsToEngine.Perform();
        //const float** inBuffers = (const float**)input;
        //float** outBuffers = (float**)output;

        int numSamples = NumSamplesPerCallback();
        int bufFrames = mWorld->mBufLength;
        int numBufs = numSamples / bufFrames;

        float* inBuses = mWorld->mAudioBus + mWorld->mNumOutputs * bufFrames;
        float* outBuses = mWorld->mAudioBus;
        int32* inTouched = mWorld->mAudioBusTouched + mWorld->mNumOutputs;
        int32* outTouched = mWorld->mAudioBusTouched;

        int bufFramePos = 0;
#ifdef SC_PA_USE_DLL
        int64 oscTime = mOSCbuftime = (uint64)((mDLL.PeriodTime() + mMaxOutputLatency) * kSecondsToOSCunits + .5);
        // 		int64 oscInc = mOSCincrement = (int64)(mOSCincrementNumerator / mDLL.SampleRate());
        int64 oscInc = mOSCincrement = (uint64)((mDLL.Period() / numBufs) * kSecondsToOSCunits + .5);
        mSmoothSampleRate = mDLL.SampleRate();
        double oscToSamples = mOSCtoSamples = mSmoothSampleRate * kOSCtoSecs /* 1/pow(2,32) */;
#else
        int64 oscTime = mOSCbuftime;
        int64 oscInc = mOSCincrement;
        double oscToSamples = mOSCtoSamples;
#endif
        // main loop
        for (int i = 0; i < numBufs; ++i, mWorld->mBufCounter++, bufFramePos += bufFrames) {
            int32 bufCounter = mWorld->mBufCounter;
            int32* tch;
/*
            // copy+touch inputs
            tch = inTouched;
            for (int k = 0; k < mInputChannelCount; ++k) {
                const float* src = inBuffers[k] + bufFramePos;
                float* dst = inBuses + k * bufFrames;
                memcpy(dst, src, bufFrames * sizeof(float));
                *tch++ = bufCounter;
            }
*/
            // run engine
            int64 schedTime;
            int64 nextTime = oscTime + oscInc;
            // DEBUG
            /*
            if (mScheduler.Ready(nextTime)) {
                double diff = (mScheduler.NextTime() - mOSCbuftime)*kOSCtoSecs;
                scprintf("rdy %.6f %.6f %.6f %.6f \n", (mScheduler.NextTime()-gStartupOSCTime) * kOSCtoSecs,
            (mOSCbuftime-gStartupOSCTime)*kOSCtoSecs, diff, (nextTime-gStartupOSCTime)*kOSCtoSecs);
            }
            */
            while ((schedTime = mScheduler.NextTime()) <= nextTime) {
                float diffTime = (float)(schedTime - oscTime) * oscToSamples + 0.5;
                float diffTimeFloor = floor(diffTime);
                world->mSampleOffset = (int)diffTimeFloor;
                world->mSubsampleOffset = diffTime - diffTimeFloor;

                if (world->mSampleOffset < 0)
                    world->mSampleOffset = 0;
                else if (world->mSampleOffset >= world->mBufLength)
                    world->mSampleOffset = world->mBufLength - 1;

                SC_ScheduledEvent event = mScheduler.Remove();
                event.Perform();
            }
            world->mSampleOffset = 0;
            world->mSubsampleOffset = 0.0f;

            World_Run(world);

            // copy touched outputs
            tch = outTouched;
            size_t pa_bytes = bufFrames * mOutputChannelCount * sizeof(float);
            float *pa_buf;
            pa_stream_begin_write(mStreamOut, (void **)&pa_buf, &pa_bytes);
            printf("CB length: %d\n", length);
            printf("pa_stream_begin_write: %d bytes\n", pa_bytes);
            for (int k = 0; k < mOutputChannelCount; ++k) {
                float* dst = (float *)pa_buf + k;// + (k * bufFrames);// + bufFramePos;
                if (tch[k] == bufCounter) {
                    const float* src = outBuses + k * bufFrames;
                    printf("dst ch %d = %p\n", k, dst);
                    //printf("src ch %d = %p\n", k, src);
                    //printf("memcpy size = %d\n", bufFrames * sizeof(float));
                    for (int n = 0; n < bufFrames; n++){
                        *dst = *src;
                        dst += mOutputChannelCount;
                        src++;
                    }
                    //memcpy(dst, src, bufFrames * sizeof(float));
                } else {
                    //memset(dst, 0, bufFrames * sizeof(float));
                    *dst = 0;
                    dst += mOutputChannelCount;
                }
            }
            pa_stream_write(mStreamOut, pa_buf, bufFrames * mOutputChannelCount * sizeof(float), NULL, 0LL, PA_SEEK_RELATIVE);

            // update buffer time
            oscTime = mOSCbuftime = nextTime;
        }
    } catch (std::exception& exc) {
        printf("SC_PulseAudioDriver: exception in real time: %s\n", exc.what());
    } catch (...) {
        printf("SC_PulseAudioDriver: unknown exception in real time\n");
    }
#if 0
    double cpuUsage = Pa_GetStreamCpuLoad(mStream) * 100.0;
    mAvgCPU = mAvgCPU + 0.1 * (cpuUsage - mAvgCPU);
    if (cpuUsage > mPeakCPU || --mPeakCounter <= 0) {
        mPeakCPU = cpuUsage;
        mPeakCounter = mMaxPeakCounter;
    }

#endif
    mAudioSync.Signal();
}

// ====================================================================
//
//


// LLP TODO: add support to choose the Output. For now, just go with the default one
void SC_PulseAudioDriver::sinkListCB(pa_context *c, const pa_sink_info *i, int eol, void *userdata)
{
    SC_PulseAudioDriver* driver = (SC_PulseAudioDriver*)userdata;
    if (eol){
        driver->sinkListReady = true;
        pa_threaded_mainloop_signal(driver->pa_ml, 0);
        return;
    }

    driver->sinkInfos[i->name] = { i->name, i->description, i->driver, i->sample_spec.channels };
}

bool SC_PulseAudioDriver::DriverSetup(int* outNumSamples, double* outSampleRate) {
    int rc;
    scprintf("In SC_PulseAudioDriver::DriverSetup\n");

    // Create a mainloop, api and context
    pa_ml = pa_threaded_mainloop_new();
    pa_mlapi = pa_threaded_mainloop_get_api(pa_ml);
    pa_ctx = pa_context_new(pa_mlapi, "SuperCollider PA");
    rc = pa_context_connect(pa_ctx, NULL, (pa_context_flags)0, NULL);
    if (rc < 0){
        PRINT_PULSEAUDIO_ERROR(pa_context_connect, rc);
        return false;
    }

    pa_context_set_state_callback(pa_ctx, SC_PulseAudioDriver::pa_state_cb, this);
    rc = pa_threaded_mainloop_start(pa_ml);    
    if (rc < 0){
        PRINT_PULSEAUDIO_ERROR(pa_threaded_mainloop_start, rc);
        return false;
    }

    pa_threaded_mainloop_lock(pa_ml);
    while (pa_ready == 0){
        pa_threaded_mainloop_wait(pa_ml);
    }
    pa_threaded_mainloop_unlock(pa_ml);
    if (pa_ready == 2) {
        PRINT_PULSEAUDIO_ERROR(pa_start, 2);
        return false;
    }

    // List the output sinks
    pa_threaded_mainloop_lock(pa_ml);
    pa_context_get_sink_info_list(pa_ctx, &SC_PulseAudioDriver::sinkListCB, this);
    while (!sinkListReady){
        pa_threaded_mainloop_wait(pa_ml);
    }
    pa_threaded_mainloop_unlock(pa_ml);

    scprintf("Available sinks: \n");
    for (const auto& sinkInfo: sinkInfos){
        scprintf("  - %s   (sink device with %d channels)\n", sinkInfo.first.c_str(), sinkInfo.second.nChannels);
    }

    // report requested devices
    //printf("\nRequested source:\n");
/*    if (mWorld->mNumInputs) {
        printf("  In (matching device %sfound):\n  - %s\n", (mDeviceInOut[0] == paNoDevice ? "NOT " : ""),
                mWorld->hw->mInDeviceName);
    }*/

    mInputChannelCount = mWorld->mNumInputs;
    mOutputChannelCount = mWorld->mNumOutputs;

    scprintf("\nRequested sink:\n");
    if (mWorld->mNumOutputs){
        if (mWorld->hw->mOutDeviceName) {
            scprintf("Name requested: %s\n", mWorld->hw->mOutDeviceName);
            auto it = sinkInfos.find(mWorld->hw->mOutDeviceName);
            if (it != sinkInfos.end()){
                scprintf("  Out (matching device found): %s\n", mWorld->hw->mOutDeviceName);
                outSinkToUse = mWorld->hw->mOutDeviceName;
                mOutputChannelCount = std::min<size_t>(mWorld->mNumOutputs, it->second.nChannels);
            }
            else{
                scprintf("  Out (matching device NOT found): %s\n", mWorld->hw->mOutDeviceName);
                mOutputChannelCount = 0;
            }
        }
        else{
            scprintf("No name requested\n");
            // select the first device
            auto it = sinkInfos.begin();
            if (it == sinkInfos.end()){
                scprintf("  No sinks found\n");
                mOutputChannelCount = 0;
            }
            else{
                scprintf("Selecting device: %s\n", it->first.c_str());
                mOutputChannelCount = std::min<size_t>(mWorld->mNumOutputs, it->second.nChannels);
                outSinkToUse = it->first;
            }
        }
    }    

/*    if (paerror != paNoError)
        PRINT_PULSEAUDIO_ERROR(Pa_OpenDefaultStream, paerror);
    return paerror == paNoError;
*/

    *outNumSamples = mWorld->mBufLength;
    if (mPreferredSampleRate)
        *outSampleRate = mPreferredSampleRate;
    else
        *outSampleRate = 44100;

    
    if (!mOutputChannelCount){
        printf("No suitable output sink found\n");
        return false;
    }

    ss.rate = *outSampleRate;
    ss.channels = mOutputChannelCount;
    ss.format = PA_SAMPLE_FLOAT32LE;
    mStreamOut = pa_stream_new(pa_ctx, "Playback", &ss, NULL);
    if (!mStreamOut) {
        printf("pa_stream_new failed\n");
        return false;
    }
    scprintf("Using %d output channels\n", mOutputChannelCount);
    pa_stream_set_write_callback(mStreamOut, &SC_PulseAudioDriver::PulseAudioCallbackStatic, this);
    //pa_stream_set_underflow_callback(playstream, stream_underflow_cb, NULL);

    return true;
}

bool SC_PulseAudioDriver::DriverStart() {
    printf("In SC_PulseAudioDriver::DriverStart\n");
    if (!mStreamOut)
        return false;
    pa_buffer_attr bufattr;

    // TODO: play with these values
    bufattr.fragsize = mWorld->mBufLength * mOutputChannelCount * sizeof(float);
    bufattr.maxlength = mWorld->mBufLength * mOutputChannelCount * sizeof(float);
    bufattr.minreq = mWorld->mBufLength * mOutputChannelCount * sizeof(float);
    bufattr.prebuf = (uint32_t)-1;
    bufattr.tlength = mWorld->mBufLength * mOutputChannelCount * sizeof(float);

    const char *outDev = nullptr;
    if (!outSinkToUse.empty()){
        outDev = outSinkToUse.c_str();
        printf("Going to use sink: %s\n", outDev);
    }

    int rc = pa_stream_connect_playback(mStreamOut, outDev, &bufattr,
                                 (pa_stream_flags)(PA_STREAM_INTERPOLATE_TIMING
                                 |PA_STREAM_ADJUST_LATENCY
                                 |PA_STREAM_AUTO_TIMING_UPDATE), NULL, NULL);
    if (rc != 0)
        PRINT_PULSEAUDIO_ERROR(pa_stream_connect_playback, rc);

    // sync times
    mPaStreamStartupTimeOSC = 0;
    mPaStreamStartupTime = 0;
    // it would be better to do the sync here, but the timeInfo in the callback is incomplete
    // mPaStreamStartupTimeOSC = GetCurrentOSCTime();
    // mPaStreamStartupTime = Pa_GetStreamTime(mStream);
#ifdef SC_PA_USE_DLL
    mDLL.Reset(mSampleRate, mNumSamplesPerCallback, SC_TIME_DLL_BW, sc_PAOSCTimeSeconds());
#endif
    return (rc == 0);
}

bool SC_PulseAudioDriver::DriverStop() {
    scprintf("In SC_PulseAudioDriver::DriverStop\n");
    if (!mStreamIn || !mStreamOut)
        return false;

    int rc = pa_stream_disconnect(mStreamIn);
    if (rcmd_af != 0)
        PRINT_PULSEAUDIO_ERROR(pa_stream_disconnect, rc);

    int rc2 = pa_stream_disconnect(mStreamOut);
    if (rc2 != 0)
        PRINT_PULSEAUDIO_ERROR(pa_stream_disconnect, rc2);

    return ((rc == 0) && (rc2 == 0));
}

void SC_PulseAudioDriver::pa_state_cb(pa_context *c, void *userdata) {
    pa_context_state_t state;
    SC_PulseAudioDriver* driver = (SC_PulseAudioDriver*)userdata;
    state = pa_context_get_state(c);
    scprintf("pa_state_cb: %d\n", (int)state);
    switch  (state) {
    // Error states
    case PA_CONTEXT_FAILED:
    case PA_CONTEXT_TERMINATED:
        driver->pa_ready = 2;
        pa_threaded_mainloop_signal(driver->pa_ml, 0);
        break;
    // Ready
    case PA_CONTEXT_READY:
        driver->pa_ready = 1;
        pa_threaded_mainloop_signal(driver->pa_ml, 0);
        break;
    // States we do not care about
    default:
        break;
    }
}