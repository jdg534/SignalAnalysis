#include "Signal.h"

#include <memory>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>


float Signal::Statistics::MeanF(const float* samples, const size_t nSamples)
{
	// prevent runtime error
	if (nSamples == 0 || samples == nullptr)
	{
		return 0.0f;
	}
	
	float sum = 0.0f;
	for (size_t i = 0; i < nSamples; ++i)
	{
		sum += samples[i];
	}
	return sum / static_cast<float>(nSamples);
}

double Signal::Statistics::MeanD(const double* samples, const size_t nSamples)
{
	// prevent runtime error
	if (nSamples == 0 || samples == nullptr)
	{
		return 0.0;
	}

	double sum = 0.0;
	for (size_t i = 0; i < nSamples; ++i)
	{
		sum += samples[i];
	}
	return sum / static_cast<double>(nSamples);
}

int16_t MeanI16(const int16_t* samples, const size_t nSamples)
{
	// prevent runtime error
	if (nSamples == 0 || samples == nullptr)
	{
		return 0;
	}

	int16_t sum = 0;
	for (size_t i = 0; i < nSamples; ++i)
	{
		sum += samples[i];
	}
	return sum / static_cast<int16_t>(nSamples);
}

float Signal::Statistics::VarienceF(const float* samples, const size_t nSamples, const float mean)
{
	// prevent runtime error, need atleast 2 samples
	if (nSamples <= 1 || samples == nullptr)
	{
		return 0.0f;
	}

	float netVarience = 0.0f;
	for (size_t i = 0; i < nSamples; ++i)
	{
		netVarience += powf(samples[i] - mean, 2);
	}
	return netVarience / static_cast<float>(nSamples - 1);
}

double Signal::Statistics::VarienceD(const double* samples, const size_t nSamples, const double mean)
{
	// prevent runtime error, need atleast 2 samples
	if (nSamples <= 1 || samples == nullptr)
	{
		return 0.0;
	}

	double netVarience = 0.0;
	for (size_t i = 0; i < nSamples; ++i)
	{
		netVarience += pow(samples[i] - mean, 2);
	}
	return netVarience / static_cast<double>(nSamples - 1);
}

float Signal::Statistics::StandardDeviationF(const float varience)
{
	return sqrtf(varience);
}

double Signal::Statistics::StandardDeviationD(const double varience)
{
	return sqrt(varience);
}

void Signal::Convolution::ConvolutionD(
	const double* inputSignalSamples, const size_t nInputSignalSamples,
	const double* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
	double* outputSignalSamples // note that outputSignalSamples is expected to be of have nInputSignalSamples + nElementInInputSignalImpulseResponse elements
	)
{
	// 0 the output buffer
	const size_t outputBufferMemSz = sizeof(double) * (nInputSignalSamples + nElementInInputSignalImpulseResponse);
	std::memset(outputSignalSamples, 0, outputBufferMemSz);

	// O(N ^ 2)
	for (size_t i = 0; i < nInputSignalSamples; ++i)
	{
		for (size_t j = 0; j < nElementInInputSignalImpulseResponse; ++j)
		{
			outputSignalSamples[i + j] = outputSignalSamples[i + j] + 
				inputSignalSamples[i] * inputSignalImpulseResponse[j];
		}
	}
}

void Signal::Convolution::ConvolutionF(
	const float* inputSignalSamples, const size_t nInputSignalSamples,
	const float* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
	float* outputSignalSamples // note that outputSignalSamples is expected to be of have nInputSignalSamples + nElementInInputSignalImpulseResponse elements
	)
{
	// 0 the output buffer
	const size_t outputBufferMemSz = sizeof(float) * (nInputSignalSamples + nElementInInputSignalImpulseResponse);
	std::memset(outputSignalSamples, 0, outputBufferMemSz);

	// O(N ^ 2)
	for (size_t i = 0; i < nInputSignalSamples; ++i)
	{
		for (size_t j = 0; j < nElementInInputSignalImpulseResponse; ++j)
		{
			outputSignalSamples[i + j] = outputSignalSamples[i + j] +
				inputSignalSamples[i] * inputSignalImpulseResponse[j];
		}
	}
}

void Signal::Convolution::RunningSumD(const double* inputSignal, const size_t nSamplesInputSignal, double* output)
{
	output[0] = inputSignal[0];
	for (size_t i = 1; i < nSamplesInputSignal; ++i)
	{
		output[i] = output[i - 1] + inputSignal[i];
	}
}

void Signal::Convolution::RunningSumF(const float* inputSignal, const size_t nSamplesInputSignal, float* output)
{
	output[0] = inputSignal[0];
	for (size_t i = 1; i < nSamplesInputSignal; ++i)
	{
		output[i] = output[i - 1] + inputSignal[i];
	}
}

void Signal::Convolution::DifferenceD(const double* inputSignal, const size_t nSamplesInInputSignal, double* output)
{
	output[0] = 0.0;
	for (size_t i = 1; i < nSamplesInInputSignal; ++i)
	{
		output[i] = inputSignal[i] - inputSignal[i - 1];
	}
}

void Signal::Convolution::DifferenceF(const float* inputSignal, const size_t nSamplesInInputSignal, float* output)
{
	output[0] = 0.0f;
	for (size_t i = 1; i < nSamplesInInputSignal; ++i)
	{
		output[i] = inputSignal[i] - inputSignal[i - 1];
	}
}

void Signal::FourierTransforms::DiscreteFourierTransformD(const double* inputSignal, const size_t inputSignalLength,
	double* outputSignalRealComponent, double* outputSignalComplexComponent)
{
	// 0 the output arrays
	const size_t outputArrayElements = (inputSignalLength / 2) + 1;
	const size_t outputBufferMemSz = sizeof(double) * outputArrayElements;
	std::memset(outputSignalRealComponent, 0, outputBufferMemSz);
	std::memset(outputSignalComplexComponent, 0, outputBufferMemSz);

	double iAsDbl = 0.0, jAsDbl;
	const double inputSignalLengthAsDbl = static_cast<double>(inputSignalLength);
	for (size_t i = 0; i < outputArrayElements; ++i, ++iAsDbl)
	{
		jAsDbl = 0.0;
		for (size_t j = 0; j < inputSignalLength; ++j, ++jAsDbl)
		{
			outputSignalRealComponent[i] += inputSignal[j] * cos(2.0 * M_PI * iAsDbl * jAsDbl / inputSignalLengthAsDbl);
			outputSignalComplexComponent[i] -= inputSignal[j] * sin(2.0 * M_PI * iAsDbl * jAsDbl / inputSignalLengthAsDbl);
		}
	}
}

void Signal::FourierTransforms::DiscreteFourierTransformMagnitudeD(double* magnitudeOutput, const double* dftRealComponent, const double* dftComplexComponent, const size_t componentArraySize)
{
	for (size_t i = 0; i < componentArraySize; ++i)
	{
		magnitudeOutput[i] = sqrt(
			(dftRealComponent[i] * dftRealComponent[i])
			+ (dftComplexComponent[i] * dftComplexComponent[i])
		);
	}
}

void Signal::FourierTransforms::InverseDiscreteFourierTransformD(double* outputSignal, const double* dftRealComponent, const double* dftComplexComponent, const size_t componentArraySize)
{
	const size_t dftComponentArraySizeInMemory = componentArraySize * sizeof(double);
	// need to convert the contents of dftRealComponent and dftComplexComponent, this requires making copies
	double* copyOfDftRealComponent = (double*) malloc(dftComponentArraySizeInMemory);
	double* copyOfDftComplexComponent = (double*)malloc(dftComponentArraySizeInMemory);
	memcpy(copyOfDftRealComponent, dftRealComponent, dftComponentArraySizeInMemory);
	memcpy(copyOfDftComplexComponent, dftComplexComponent, dftComponentArraySizeInMemory);

	const size_t outSigArraySize = (componentArraySize - 1) * 2;

	// convert the component array's samples be: sample[i] /= (outSigArraySize / 2)
	const double outSigArraySizeOverTwo = static_cast<double>(outSigArraySize) / 2.0;
	for (size_t i = 0; i < componentArraySize; ++i)
	{
		copyOfDftRealComponent[i] /= outSigArraySizeOverTwo;
		copyOfDftComplexComponent[i] /= outSigArraySizeOverTwo;
	}

	// 0 the output array	
	const size_t outSizeMemSize = outSigArraySize * sizeof(double);
	std::memset(outputSignal, 0, outSizeMemSize);

	const double iDftLength = static_cast<double>(outSigArraySize);
	double dblDftSampleIndex = 0.0, dblOutputSigSampleIndex = 0.0;
	for (size_t dftSampleIndex = 0; dftSampleIndex < componentArraySize; ++dftSampleIndex, ++dblDftSampleIndex)
	{
		dblOutputSigSampleIndex = 0.0;
		for (size_t outputSignalSampleIndex = 0; outputSignalSampleIndex < outSigArraySize; ++outputSignalSampleIndex, ++dblOutputSigSampleIndex)
		{
			outputSignal[outputSignalSampleIndex] += copyOfDftRealComponent[dftSampleIndex] * cos(2.0 * M_PI * dblDftSampleIndex * dblOutputSigSampleIndex / iDftLength);
			outputSignal[outputSignalSampleIndex] += copyOfDftComplexComponent[dftSampleIndex] * sin(2.0 * M_PI * dblDftSampleIndex * dblOutputSigSampleIndex / iDftLength);
		}
	}

	// clean up don't want to leak memory
	free(copyOfDftRealComponent);
	free(copyOfDftComplexComponent);
}

void Signal::Filters::MovingAverageSubsquentPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage)
{
	const double nPointsToAvgAsDbl = static_cast<double>(nPointsToAverage);
	double average;

	for (size_t i = 0; i < inputSignalLength; ++i)
	{
		average = 0.0;
		for (size_t j = i; j < inputSignalLength && j < nPointsToAverage + i; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsDbl;
		output[i] = average;
	}
}

void Signal::Filters::MovingAverageSymetricallyChosenPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage)
{
	const double nPointsToAvgAsDbl = static_cast<double>(nPointsToAverage);
	double average;

	const size_t halfNumPointsToAverage = nPointsToAverage / 2;
	const int inputSigLenAsSignedInt = static_cast<int>(inputSignalLength);

	for (int i = 0; i < inputSignalLength; ++i)
	{
		average = 0.0;
		const int posibleStartingIndex = i - halfNumPointsToAverage;
		const int sigAvgStartingPoint = posibleStartingIndex >= 0 ? posibleStartingIndex : 0; // start at i - halfNumPointsToAverage or the start of the array if we can't go back that far
		for (size_t j = sigAvgStartingPoint; j < inputSignalLength && j < sigAvgStartingPoint + nPointsToAverage; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsDbl;
		output[i] = average;
	}
}

void Signal::Filters::Windowed::HammingD(double* outputWindow, const size_t inputSignalLength)
{
	// equation = w[i] = 0.54 - 0.46 * cos(2*PI*i / M)
	// i = sample index
	// M = number of samples

	const double dblInputSigLength = static_cast<double>(inputSignalLength);

	double iterAsDbl = 0.0;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsDbl)
	{
		outputWindow[i] = 0.54 - 0.46 * cos(2.0 * M_PI * iterAsDbl / dblInputSigLength);
	}
}

void Signal::Filters::Windowed::BlackmanD(double* outputWindow, const size_t inputSignalLength)
{
	// do Blackman window next
	// equation = w[i] = 0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M)
	// i = sample index
	// M = number of samples

	const double dblInputSigLength = static_cast<double>(inputSignalLength);

	double iterAsDbl = 0.0;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsDbl)
	{
		outputWindow[i] = 0.42 - 0.5 * cos(2.0 * M_PI * iterAsDbl / dblInputSigLength) + 0.08 * cos(4 * M_PI * iterAsDbl / dblInputSigLength);
	}
}

void Signal::Filters::Windowed::SyncLowPassD(double* filterOutput, const size_t filterOutputSize, double cutoffFrequency)
{
	const int filterOutputSizeI = static_cast<int>(filterOutputSize); // need signed values
	const int halfFilterLength = filterOutputSizeI / 2;

	double iDbl = 0.0;
	const double filterOutputSizeDbl = static_cast<double>(filterOutputSize);

	for (int i = 0; i < filterOutputSizeI; ++i, ++iDbl)
	{
		if (i - halfFilterLength == 0) // half way point is a special value
		{
			filterOutput[i] = 2.0 * M_PI * cutoffFrequency;
		}
		else
		{
			filterOutput[i] = sin(2.0 * M_PI * cutoffFrequency * (static_cast<double>(i - halfFilterLength)))
				/ (static_cast<double>(i - halfFilterLength));
			filterOutput[i] *= (0.54 - 0.46 * cos(2.0 * M_PI * iDbl / filterOutputSizeDbl));
		}
	}
}