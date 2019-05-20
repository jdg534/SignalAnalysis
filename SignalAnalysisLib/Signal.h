#ifndef _SIGNAL_H_
#define _SIGNAL_H_

// will go for name spaces as was the case in the Josh Maths library
namespace Signal
{
	namespace Statistics
	{
		float MeanF(const float* samples, const size_t nSamples);
		double MeanD(const double* samples, const size_t nSamples);
	}
}

#endif // _SIGNAL_H_