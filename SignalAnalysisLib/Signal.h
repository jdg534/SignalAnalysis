#ifndef _SIGNAL_H_
#define _SIGNAL_H_

// will go for name spaces as was the case in the Josh Maths library
namespace Signal
{
	namespace Statistics
	{
		float MeanF(const float* samples, const size_t nSamples);
		double MeanD(const double* samples, const size_t nSamples);
		float VarienceF(const float* samples, const size_t nSamples, const float mean);
		double VarienceD(const double* samples, const size_t nSamples, const double mean);
		float StandardDeviationF(const float varience);
		double StandardDeviationD(const double varience);
	}
	namespace Convolution
	{
		void ConvolutionD(const double* inputSignalSamples, const size_t nInputSignalSamples,
			double* outputSignalSamples, // note that outputSignalSamples is assumed to be the same size as nInputSignalSamples
			const double* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse);
	}
}

#endif // _SIGNAL_H_