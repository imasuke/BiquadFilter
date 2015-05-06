#ifndef _SIGNAL_CLASS_TEMPLATE_HEADER_
#define _SIGNAL_CLASS_TEMPLATE_HEADER_


#include <iostream>
#include <vector>

template<class T> class SignalT
{
public:
    // sample value and time
    typedef struct SampleWithTime_ {
        double time_;
        T sample_value_;
    } SampleWithTime;
    typedef std::vector<SampleWithTime> Samples;
    // spectrum
    typedef struct Complex_ {
        double real;
        double imag;
        Complex_(double real_, double imag_)
            : real(real_), imag(imag_){ }
    } Complex;
    typedef std::vector<Complex> Spectrum;

protected:
    std::size_t spFreq_;    // sampling frequency
    std::size_t qBit_;      // quantized bits
    std::size_t length_;    // signal length
    Samples     samples_;   // raw sample data

public:
    SignalT(void)
        : spFreq_(0), qBit_(0), length_(0)
    {
        samples_.clear();
    }
    SignalT(std::size_t spFreq, std::size_t qBit, std::size_t length)
        : spFreq_(spFreq), qBit_(qBit), length_(length)
    {
        samples_.resize(length_);
    }
    SignalT(const SignalT<T> &signal)
    {
        spFreq_  = signal.spFreq_;
        qBit_    = signal.qBit_;
        length_  = signal.length_;
        samples_ = signal.samples_;
    }
    SignalT<T> &operator=(const SignalT<T> &signal)
    {
        spFreq_  = signal.spFreq_;
        qBit_    = signal.qBit_;
        length_  = signal.length_;
        samples_ = signal.samples_;
        return (*this);
    }
    ~SignalT(void)
    {

    }

    void allocate(std::size_t spFreq, std::size_t qBit, std::size_t length)
    {
        spFreq_ = spFreq;
        qBit_   = qBit;
        length_ = length;
        samples_.resize(length);
    }

    T operator[](std::size_t n) const
    {
        return samples_[n].sample_value_;
    }
    T &operator[](std::size_t n)
    {
        return samples_[n].sample_value_;
    }
    
    double dt(void)
    {
        return 1.0 / static_cast<double>(spFreq_);
    }
    double df(void)
    {
        return spFreq_ / static_cast<double>(samples_.size());
    }

    std::size_t spFreq(void) const
    {
      return spFreq_;
    }
    std::size_t qBit(void) const 
    {
      return qBit_;
    }
    std::size_t length(void) const
    {
      return length_;
    }
    Samples samples(void) const
    {
      return samples_;
    }

};

#endif // _SIGNAL_TEMPLATE_HEADER_