
//  $Id: front_end.C,v 1.2 2016/01/23 03:15:23 stanchen Exp $


#include "front_end.H"
#include "math.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** Module for doing windowing. **/
void FrontEnd::do_window(const matrix<double>& inFeats,
    matrix<double>& outFeats) const
    {
    //  Get parameters.
    //  Input samples per second.
    double sampleRate = get_float_param(m_params, "window.sample_rate",
        20000.0);
    //  Output frames per second.
    double framesPerSec = get_float_param(m_params, "window.frames_per_sec",
        100.0);
    //  Width of each window, in seconds.
    double windowWidth = get_float_param(m_params, "window.window_size",
        0.025);
    //  Whether to do Hamming or rectangular windowing.
    bool doHamming = get_bool_param(m_params, "window.hamming", true);

    //  Get number of input samples.
    int inSampCnt = inFeats.size1();
    if (inFeats.size2() != 1)
        throw runtime_error("Windowing expected vector input.");

    //  Input sampling period in seconds.
    double samplePeriod = 1.0 / sampleRate;
    //  Output frame period, in seconds.
    double framePeriod = 1.0 / framesPerSec;
    //  Number of samples per window.
    int sampPerWindow = (int) (windowWidth / samplePeriod + 0.5);
    //  Number of samples to shift between each window.
    int sampShift = (int) (framePeriod / samplePeriod + 0.5);
    //  Number of output frames.
    int outFrameCnt = (inSampCnt - sampPerWindow) / sampShift + 1;

    //  Allocate output matrix and fill with zeros.
    outFeats.resize(outFrameCnt, sampPerWindow);
    outFeats.clear();
    
    int window_left;
    // windowing
    for (int i=0;i<outFrameCnt;i++){
        window_left=i*sampShift;
        for (int j=window_left;j<window_left+sampPerWindow;j++){
            if (doHamming==false){
                outFeats(i,j-window_left)=inFeats(j,0);
            }
            else{
                outFeats(i,j-window_left)=(0.54-0.46*cos(2*M_PI*(j-window_left)/(sampPerWindow-1)))*inFeats(j,0);
            }
        }
    }

    //  BEGIN_LAB
    //
    //  Input:
    //      "inFeats", a matrix containing a single column holding the
    //      input samples for an utterance.  Each row holds a single sample.
    //
    //      inFeats(0 .. (inSampCnt - 1), 0)
    //
    //  Output:
    //      "outFeats", which should contain the result of windowing.
    //
    //      outFeats(0 .. (outFrameCnt - 1), 0 .. (sampPerWindow - 1))
    //
    //      Each row corresponds to a frame, and should hold the
    //      windowed samples for that frame.
    //      It has already been allocated to be of the correct size.
    //      If the boolean "doHamming" is true, then a Hamming
    //      window should be applied; otherwise, a rectangular
    //      window should be used.
    //
    //  See "inSampCnt", "sampPerWindow", "sampShift", and "outFrameCnt"
    //  above for quantities you may (or may not) need for this computation.
    //
    //  When accessing matrices such as "inFeats" and "outFeats",
    //  use a syntax like "inFeats(frmIdx, dimIdx)" to access elements;
    //  using square brackets as in normal C arrays won't work.
    




    //  END_LAB
    }

/** Module for doing FFT. **/
void FrontEnd::do_fft(const matrix<double>& inFeats,
    matrix<double>& outFeats) const
    {
    //  Make output dimension the smallest power of 2 at least as
    //  large as input dimension.
    int inFrameCnt = inFeats.size1();
    int inDimCnt = inFeats.size2();
    int outDimCnt = 2;
    while (outDimCnt < inDimCnt)
        outDimCnt *= 2;

    //  Allocate output matrix and fill with zeros.
    outFeats.resize(inFrameCnt, outDimCnt);
    outFeats.clear();

    //  Input:
    //      "inFeats", a matrix with each row holding the windowed
    //      values for that frame.
    //
    //      inFeats(0 .. (inFrameCnt - 1), 0 .. (inDimCnt - 1))
    //
    //  Output:
    //      "outFeats", where an FFT should be applied to each
    //      row/frame of "inFeats".  
    //
    //      outFeats(0 .. (inFrameCnt - 1), 0 .. (outDimCnt - 1))
    //
    //      For a given row/frame "frmIdx", the real and imaginary
    //      parts of the FFT value for frequency i/(outDimCnt*T)
    //      where T is the sample period are held in
    //      outFeats(frmIdx, 2*i) and outFeats(frmIdx, 2*i+1),
    //      respectively.

    vector<double> fftBuf;
    for (int frmIdx = 0; frmIdx < inFrameCnt; ++frmIdx)
        {
        copy_matrix_row_to_vector(inFeats, frmIdx, fftBuf);
        //  Pad window with zeros, if needed.
        fftBuf.resize(outDimCnt, 0.0);
        real_fft(fftBuf);
        copy_vector_to_matrix_row(fftBuf, outFeats, frmIdx);
        }
    }

/** Module for mel binning. **/
void FrontEnd::do_melbin(const matrix<double>& inFeats,
    matrix<double>& outFeats) const
    {
    //  Number of mel bins to make.
    int numBins = get_int_param(m_params, "melbin.bins", 26);
    //  Whether to take log of output or not.
    bool doLog = get_bool_param(m_params, "melbin.log", true);
    //  Input samples per second.
    double sampleRate = get_float_param(m_params, "window.sample_rate",
        20000.0);
    double samplePeriod = 1.0 / sampleRate;

    //  Retrieve number of frames and dimension of input feature vectors.
    int inFrameCnt = inFeats.size1();
    int inDimCnt = inFeats.size2();
    int outDimCnt = numBins;

    //  Allocate output matrix and fill with zeros.
    outFeats.resize(inFrameCnt, outDimCnt);
    outFeats.clear();
    
    double X_f;
    double freq;
    double H;

    for (int i=0;i<inFrameCnt;i++){
        for (int j=0;j<numBins;j++){
            for (int k=0;k<inDimCnt;k+=2){
                X_f=sqrt(pow(inFeats(i,k),2)+pow(inFeats(i,k+1),2));
                freq=mel_func(k/(inDimCnt*samplePeriod)/2.0);
                H=0.0;
                if (freq<=func(j+1) && freq>=func(j)){
                    H=(freq-func(j))/(func(j+1)-func(j))/1.0;
                }
                if (freq<=func(j+2) && freq>=func(j+1)){
                    H=(func(j+2)-freq)/(func(j+2)-func(j+1))/1.0;
                }
                outFeats(i,j)+=X_f*H;
            }

        }

    }

    if (doLog==1) {
        for (int i=0;i<inFrameCnt;i++){
            for (int j=0;j<numBins;j++){
                if (log(outFeats(i,j))>0){
                    outFeats(i,j)=log(outFeats(i,j));
                }
                else {
                    outFeats(i,j)=0;
                }
            }
        }
    }
    //  BEGIN_LAB
    //
    //  Input:
    //      "inFeats", holding the output of a real FFT.
    //
    //      inFeats(0 .. (inFrameCnt - 1), 0 .. (inDimCnt - 1))
    //
    //  Output:
    //      "outFeats", which should contain the result of
    //      mel-binning.
    //
    //      outFeats(0 .. (inFrameCnt - 1), 0 .. (outDimCnt - 1))
    //
    //      If the boolean "doLog" is true,
    //      then each value should be replaced with its natural
    //      logarithm, or 0 if its logarithm is negative.
    //      "outFeats" has been allocated to be of the correct size.
    //
    //  See "inFrameCnt", "inDimCnt", "outDimCnt", and "samplePeriod"
    //  above for quantities you will need for this computation.
    //  END_LAB
    }
double FrontEnd::func(int m) const
{
    int numBins = get_int_param(m_params, "melbin.bins", 26);
    double sampleRate = get_float_param(m_params, "window.sample_rate",
        20000.0);
    double samplePeriod = 1.0 / sampleRate;
    return mel_func(1.0/(2*samplePeriod))/(numBins+1)*m;
}
double FrontEnd::mel_func(double f) const
{
    return 1127*log(1.0+f/700.0);
}

/** Module for doing discrete cosine transform. **/
void FrontEnd::do_dct(const matrix<double>& inFeats,
    matrix<double>& outFeats) const
    {
    //  Number of DCT coefficients to output.
    int numCoeffs = get_int_param(m_params, "dct.coeffs", 12);
    int inFrameCnt = inFeats.size1();
    int inDimCnt = inFeats.size2();
    int outDimCnt = numCoeffs;

    //  Allocate output matrix and fill with zeros.
    outFeats.resize(inFrameCnt, outDimCnt);
    outFeats.clear();
    
    for (int i=0;i<inFrameCnt;i++){
        for (int j=0;j<outDimCnt;j++){
           for (int k=0;k<inDimCnt;k++){
            outFeats(i,j)+=sqrt(2.0/inDimCnt)*inFeats(i,k)*cos(M_PI*(j+1)/inDimCnt*(k+0.5));
           }
        }
    }
    //  BEGIN_LAB
    //
    //  Input:
    //      The matrix "inFeats", holding the output of mel-binning.
    //
    //      inFeats(0 .. (inFrameCnt - 1), 0 .. (inDimCnt - 1))
    //
    //  Output:
    //      The matrix "outFeats", which should contain the result of
    //      applying the DCT.
    //
    //      outFeats(0 .. (inFrameCnt - 1), 0 .. (outDimCnt - 1))
    //
    //      "outFeats" has been allocated to be of the correct size.
    //
    //  See "inFrameCnt", "inDimCnt", and "outDimCnt" above
    //  for quantities you will need for this computation.





    //  END_LAB
    }

/** Main signal processing routine.
*   Calls each signal processing module in turn, unless
*   parameter says not to.
**/
void FrontEnd::get_feats(const matrix<double>& inAudio,
    matrix<double>& outFeats) const
    {
    if (get_bool_param(m_params, "frontend.null", false))
        {
        outFeats = inAudio;
        return;
        }
    matrix<double> curFeats(inAudio);
    if (get_bool_param(m_params, "frontend.window", true))
        {
        do_window(curFeats, outFeats);
        outFeats.swap(curFeats);
        }
    if (get_bool_param(m_params, "frontend.fft", true))
        {
        do_fft(curFeats, outFeats);
        outFeats.swap(curFeats);
        }
    if (get_bool_param(m_params, "frontend.melbin", true))
        {
        do_melbin(curFeats, outFeats);
        outFeats.swap(curFeats);
        }
    if (get_bool_param(m_params, "frontend.dct", true))
        {
        do_dct(curFeats, outFeats);
        outFeats.swap(curFeats);
        }
    outFeats.swap(curFeats);
    }


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


