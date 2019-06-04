#ifndef _SIGNAL_H_
#define _SIGNAL_H_

#include <cstdint>

// will go for name spaces as was the case in the Josh Maths library
namespace Signal
{
	// support float, double, S16

	namespace Statistics
	{
		float MeanF(const float* samples, const size_t nSamples);
		double MeanD(const double* samples, const size_t nSamples);
		int16_t MeanI16(const int16_t* samples, const size_t nSamples);
		float VarienceF(const float* samples, const size_t nSamples, const float mean);
		double VarienceD(const double* samples, const size_t nSamples, const double mean);
		int16_t VarienceI16(const int16_t * samples, const size_t nSamples, const int16_t  mean);
		float StandardDeviationF(const float varience);
		double StandardDeviationD(const double varience);
		int16_t StandardDeviationI16(const int16_t varience);
	}
	namespace Convolution
	{
		void ConvolutionD(
			const double* inputSignalSamples, const size_t nInputSignalSamples,
			const double* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
			double* outputSignalSamples // note that outputSignalSamples is expected to have (nInputSignalSamples + nElementInInputSignalImpulseResponse) elements
			);
		void ConvolutionF(
			const float* inputSignalSamples, const size_t nInputSignalSamples,
			const float* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
			float* outputSignalSamples // note that outputSignalSamples is expected to have (nInputSignalSamples + nElementInInputSignalImpulseResponse) elements
			);
		void ConvolutionI16(
			const int16_t* inputSignalSamples, const size_t nInputSignalSamples,
			const int16_t* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
			int16_t* outputSignalSamples // note that outputSignalSamples is expected to have (nInputSignalSamples + nElementInInputSignalImpulseResponse) elements
		);

		// output is expected to be nSamplesInputSignal in size
		void RunningSumD(const double* inputSignal, const size_t nSamplesInputSignal, double* output);
		void RunningSumF(const float* inputSignal, const size_t nSamplesInputSignal, float* output);
		void RunningSumI16(const int16_t* inputSignal, const size_t nSamplesInputSignal, int16_t* output);

		// output expected to be the same size as inputSignal
		void DifferenceD(const double* inputSignal, const size_t nSamplesInInputSignal, double* output); 
		void DifferenceF(const float* inputSignal, const size_t nSamplesInInputSignal, float* output);
		void DifferenceI16(const int16_t* inputSignal, const size_t nSamplesInInputSignal, int16_t* output);
	}

	namespace FourierTransforms
	{
		void DiscreteFourierTransformD(const double* inputSignal, const size_t inputSignalLength, 
			double* outputSignalRealComponent,
			double* outputSignalComplexComponent);

		void DiscreteFourierTransformMagnitudeD(double* magnitudeOutput, const double* dftRealComponent, const double* dftComplexComponent, const size_t componentArraySize);

		void InverseDiscreteFourierTransformD(double* outputSignal, const double* dftRealComponent, const double* dftComplexComponent, const size_t componentArraySize);
	}

	namespace Filters
	{
		void MovingAverageSubsquentPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage);
		void MovingAverageSymetricallyChosenPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage);

		// note the recursive is the DSP definition for the term (RecursiveMovingAverageD doesn't call itself)
		void RecursiveMovingAverageD(const double* intputSignal, const size_t inputSignalLength, double* outputSignal, const size_t nPointsToAverage);

		namespace Windowed
		{
			void HammingD(double* outputWindow, const size_t inputSignalLength);
			void BlackmanD(double* outputWindow, const size_t inputSignalLength);
			void SyncLowPassD(double* filterOutput, const size_t filterOutputSize, double cutoffFrequency);
		}
	}
}

#endif // _SIGNAL_H_