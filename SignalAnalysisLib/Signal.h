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
		void ConvolutionD(
			const double* inputSignalSamples, const size_t nInputSignalSamples,
			const double* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
			double* outputSignalSamples // note that outputSignalSamples is expected to be of have nInputSignalSamples + nElementInInputSignalImpulseResponse elements
			);
		void ConvolutionF(
			const float* inputSignalSamples, const size_t nInputSignalSamples,
			const float* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
			float* outputSignalSamples // note that outputSignalSamples is expected to be of have nInputSignalSamples + nElementInInputSignalImpulseResponse elements
			);

		// output is expected to be nSamplesInputSignal in size
		void RunningSumD(const double* inputSignal, const size_t nSamplesInputSignal, double* output);
		void RunningSumF(const float* inputSignal, const size_t nSamplesInputSignal, float* output);

		// difference
		void DifferenceD(const double* inputSignal, const size_t nSamplesInInputSignal, double* output); // output expected to be the same size as inputSignal
		void DifferenceF(const float* inputSignal, const size_t nSamplesInInputSignal, float* output); // output expected to be the same size as inputSignal
	}

	namespace FourierTransforms
	{
		void DiscreteFourierTransformD(const double* inputSignal, const size_t inputSignalLength, 
			double* outputSignalRealComponent,
			double* outputSignalComplexComponent);

		// note magnitudeOutput is expected to be (componentArraySize / 2) in length
		void DiscreteFourierTransformMagnitudeD(double* magnitudeOutput, const double* dftRealComponent, const double* dftComplexComponent, const size_t componentArraySize);

		void InverseDiscreteFourierTransformD(double* outputSignal, const double* dftRealComponent, const double* dftComplexComponent, const size_t componentArraySize);
	}

	namespace Filters
	{
		void MovingAverageSubsquentPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage);
		void MovingAverageSymetricallyChosenPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage);

		namespace Windowed
		{
			void HammingD(const double* inputSignal, const size_t inputSignalLength, double* outputWindow);
			void BlackmanD(const double* inputSignal, const size_t inputSignalLength, double* outputWindow);
			void SyncLowPassD(double* filterOutput, const size_t filterOutputSize, double cutoffFrequency);
		}
	}
}

#endif // _SIGNAL_H_