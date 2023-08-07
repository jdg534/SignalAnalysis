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

int16_t Signal::Statistics::MeanI16(const int16_t* samples, const size_t nSamples)
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
	// prevent runtime error, need at least 2 samples
	if (nSamples <= 1 || samples == nullptr)
	{
		return 0.0f;
	}

	float netVarience = 0.0f;
	float deltaFromMean = 0.0f;
	for (size_t i = 0; i < nSamples; ++i)
	{
		deltaFromMean = samples[i] - mean;
		netVarience += (deltaFromMean * deltaFromMean);
	}
	return netVarience / static_cast<float>(nSamples - 1);
}

double Signal::Statistics::VarienceD(const double* samples, const size_t nSamples, const double mean)
{
	// prevent runtime error, need at least 2 samples
	if (nSamples <= 1 || samples == nullptr)
	{
		return 0.0;
	}

	double netVarience = 0.0;
	double deltaFromMean = 0.0;
	for (size_t i = 0; i < nSamples; ++i)
	{
		deltaFromMean = samples[i] - mean;
		netVarience += (deltaFromMean * deltaFromMean);
	}
	return netVarience / static_cast<double>(nSamples - 1);
}

int16_t Signal::Statistics::VarienceI16(const int16_t * samples, const size_t nSamples, const int16_t  mean)
{
	// prevent runtime error, need at least 2 samples
	if (nSamples <= 1 || samples == nullptr)
	{
		return 0;
	}

	int16_t netVarience = 0;
	int16_t deltaFromMean = 0;
	for (size_t i = 0; i < nSamples; ++i)
	{
		deltaFromMean = samples[i] - mean;
		netVarience += (deltaFromMean * deltaFromMean);
	}
	return netVarience / static_cast<int16_t>(nSamples - 1);
}

float Signal::Statistics::StandardDeviationF(const float varience)
{
	return sqrtf(varience);
}

double Signal::Statistics::StandardDeviationD(const double varience)
{
	return sqrt(varience);
}

int16_t Signal::Statistics::StandardDeviationI16(const int16_t varience)
{
	return static_cast<int16_t>(sqrt(varience));
}

void Signal::Convolution::ConvolutionD(
	const double* inputSignalSamples, const size_t nInputSignalSamples,
	const double* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
	double* outputSignalSamples // note that outputSignalSamples is expected to have (nInputSignalSamples + nElementInInputSignalImpulseResponse) elements
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
	float* outputSignalSamples // note that outputSignalSamples is expected to have (nInputSignalSamples + nElementInInputSignalImpulseResponse) elements
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

void Signal::Convolution::ConvolutionI16(
	const int16_t* inputSignalSamples, const size_t nInputSignalSamples,
	const int16_t* inputSignalImpulseResponse, const size_t nElementInInputSignalImpulseResponse,
	int16_t* outputSignalSamples // note that outputSignalSamples is expected to have (nInputSignalSamples + nElementInInputSignalImpulseResponse) elements
)
{
	// 0 the output buffer
	const size_t outputBufferMemSz = sizeof(int16_t) * (nInputSignalSamples + nElementInInputSignalImpulseResponse);
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

void Signal::Convolution::RunningSumI16(const int16_t* inputSignal, const size_t nSamplesInputSignal, int16_t* output)
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

void Signal::Convolution::DifferenceI16(const int16_t* inputSignal, const size_t nSamplesInInputSignal, int16_t* output)
{
	output[0] = 0;
	for (size_t i = 1; i < nSamplesInInputSignal; ++i)
	{
		output[i] = inputSignal[i] - inputSignal[i - 1];
	}
}

size_t Signal::FourierTransforms::OutputSamplesLength(const size_t inputSignalLength)
{
	return (inputSignalLength / 2) + 1;
}

void Signal::FourierTransforms::DiscreteFourierTransformD(const double* inputSignal, const size_t inputSignalLength,
	double* outputSignalRealComponent, double* outputSignalComplexComponent)
{
	// 0 the output arrays
	const size_t outputArrayElements = OutputSamplesLength(inputSignalLength);
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

void Signal::FourierTransforms::DiscreteFourierTransformF(const float* inputSignal, const size_t inputSignalLength,
	float* outputSignalRealComponent, float* outputSignalComplexComponent)
{
	// 0 the output arrays
	const size_t outputArrayElements = OutputSamplesLength(inputSignalLength);
	const size_t outputBufferMemSz = sizeof(float) * outputArrayElements;
	std::memset(outputSignalRealComponent, 0, outputBufferMemSz);
	std::memset(outputSignalComplexComponent, 0, outputBufferMemSz);

	float iAsFlt = 0.0f, jAsFlt;
	const float inputSignalLengthAsFlt = static_cast<float>(inputSignalLength);
	const float fltPi = static_cast<float>(M_PI);
	for (size_t i = 0; i < outputArrayElements; ++i, ++iAsFlt)
	{
		jAsFlt = 0.0f;
		for (size_t j = 0; j < inputSignalLength; ++j, ++jAsFlt)
		{
			outputSignalRealComponent[i] += inputSignal[j] * cosf(2.0f * fltPi * iAsFlt * jAsFlt / inputSignalLengthAsFlt);
			outputSignalComplexComponent[i] -= inputSignal[j] * sinf(2.0f * fltPi * iAsFlt * jAsFlt / inputSignalLengthAsFlt);
		}
	}
}

void Signal::FourierTransforms::DiscreteFourierTransformI16(const int16_t* inputSignal, const size_t inputSignalLength,
	int16_t* outputSignalRealComponent, int16_t* outputSignalComplexComponent)
{
	// 0 the output arrays
	const size_t outputArrayElements = OutputSamplesLength(inputSignalLength);
	const size_t outputBufferMemSz = sizeof(int16_t) * outputArrayElements;
	std::memset(outputSignalRealComponent, 0, outputBufferMemSz);
	std::memset(outputSignalComplexComponent, 0, outputBufferMemSz);

	float iAsFlt = 0.0f, jAsFlt;
	const float inputSignalLengthAsFlt = static_cast<float>(inputSignalLength);
	const float fltPi = static_cast<float>(M_PI);

	// will need to put the int16 range values into float range values
	// float range (logical) = -1.0 to 1.0 (no capping)
	// INT16 range to -32768 t0 +32767
	constexpr float scaleToFloatRange = 1.0f / 32767.0f;
	constexpr float scaleToInt16Range = 32767.0f / 1.0f;

	float fltRealComponent = 0.0f;
	float fltComplexComponent = 0.0f;

	for (size_t i = 0; i < outputArrayElements; ++i, ++iAsFlt)
	{
		jAsFlt = 0.0f;
		for (size_t j = 0; j < inputSignalLength; ++j, ++jAsFlt)
		{
			fltRealComponent = static_cast<float>(inputSignal[j]);
			fltComplexComponent = static_cast<float>(inputSignal[j]);
			// put value into float range
			fltRealComponent *= scaleToFloatRange;
			fltComplexComponent *= scaleToFloatRange;
			// calculate with float range logic
			fltRealComponent *= cosf(2.0f * fltPi * iAsFlt * jAsFlt / inputSignalLengthAsFlt);
			fltComplexComponent *= sinf(2.0f * fltPi * iAsFlt * jAsFlt / inputSignalLengthAsFlt);
			// put calculated value into INT16 range
			fltRealComponent *= scaleToInt16Range;
			fltComplexComponent *= scaleToInt16Range;
			outputSignalRealComponent[i] += static_cast<int16_t>(fltRealComponent);
			outputSignalComplexComponent[i] -= static_cast<int16_t>(fltComplexComponent);
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

void Signal::FourierTransforms::DiscreteFourierTransformMagnitudeF(float* magnitudeOutput, const float* dftRealComponent, const float* dftComplexComponent, const size_t componentArraySize)
{
	for (size_t i = 0; i < componentArraySize; ++i)
	{
		magnitudeOutput[i] = sqrtf(
			(dftRealComponent[i] * dftRealComponent[i])
			+ (dftComplexComponent[i] * dftComplexComponent[i])
		);
	}
}

void Signal::FourierTransforms::DiscreteFourierTransformMagnitudei16(int16_t* magnitudeOutput, const int16_t* dftRealComponent, const int16_t* dftComplexComponent, const size_t componentArraySize)
{
	int16_t squaredComps;
	for (size_t i = 0; i < componentArraySize; ++i)
	{
		squaredComps = (dftRealComponent[i] * dftRealComponent[i])
			+ (dftComplexComponent[i] * dftComplexComponent[i]);
		magnitudeOutput[i] = static_cast<int16_t>(sqrt(squaredComps));
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

void Signal::FourierTransforms::InverseDiscreteFourierTransformF(float* outputSignal, const float* dftRealComponent, const float* dftComplexComponent, const size_t componentArraySize)
{
	const size_t dftComponentArraySizeInMemory = componentArraySize * sizeof(float);
	// need to convert the contents of dftRealComponent and dftComplexComponent, this requires making copies
	float* copyOfDftRealComponent = (float*)malloc(dftComponentArraySizeInMemory);
	float* copyOfDftComplexComponent = (float*)malloc(dftComponentArraySizeInMemory);
	memcpy(copyOfDftRealComponent, dftRealComponent, dftComponentArraySizeInMemory);
	memcpy(copyOfDftComplexComponent, dftComplexComponent, dftComponentArraySizeInMemory);

	const size_t outSigArraySize = (componentArraySize - 1) * 2;

	// convert the component array's samples be: sample[i] /= (outSigArraySize / 2)
	const float outSigArraySizeOverTwo = static_cast<float>(outSigArraySize) / 2.0f;
	for (size_t i = 0; i < componentArraySize; ++i)
	{
		copyOfDftRealComponent[i] /= outSigArraySizeOverTwo;
		copyOfDftComplexComponent[i] /= outSigArraySizeOverTwo;
	}

	// 0 the output array	
	const size_t outSizeMemSize = outSigArraySize * sizeof(float);
	std::memset(outputSignal, 0, outSizeMemSize);

	const float iDftLength = static_cast<float>(outSigArraySize);
	float dblDftSampleIndex = 0.0f, dblOutputSigSampleIndex = 0.0f;
	const float fltPi = static_cast<float>(M_PI);
	for (size_t dftSampleIndex = 0; dftSampleIndex < componentArraySize; ++dftSampleIndex, ++dblDftSampleIndex)
	{
		dblOutputSigSampleIndex = 0.0f;
		for (size_t outputSignalSampleIndex = 0; outputSignalSampleIndex < outSigArraySize; ++outputSignalSampleIndex, ++dblOutputSigSampleIndex)
		{
			outputSignal[outputSignalSampleIndex] += copyOfDftRealComponent[dftSampleIndex] * cos(2.0f * fltPi * dblDftSampleIndex * dblOutputSigSampleIndex / iDftLength);
			outputSignal[outputSignalSampleIndex] += copyOfDftComplexComponent[dftSampleIndex] * sin(2.0f * fltPi * dblDftSampleIndex * dblOutputSigSampleIndex / iDftLength);
		}
	}

	// clean up don't want to leak memory
	free(copyOfDftRealComponent);
	free(copyOfDftComplexComponent);
}

void Signal::FourierTransforms::InverseDiscreteFourierTransformI16(int16_t* outputSignal, const int16_t* dftRealComponent, const int16_t* dftComplexComponent, const size_t componentArraySize)
{
	// just to the double logic and convert numbers to the correct value range

	const size_t dftComponentArraySizeInMemory = componentArraySize * sizeof(int16_t);
	// need to convert the contents of dftRealComponent and dftComplexComponent, this requires making copies
	int16_t* copyOfDftRealComponent = (int16_t*)malloc(dftComponentArraySizeInMemory);
	int16_t* copyOfDftComplexComponent = (int16_t*)malloc(dftComponentArraySizeInMemory);
	memcpy(copyOfDftRealComponent, dftRealComponent, dftComponentArraySizeInMemory);
	memcpy(copyOfDftComplexComponent, dftComplexComponent, dftComponentArraySizeInMemory);

	const size_t outSigArraySize = (componentArraySize - 1) * 2;

	// convert the component array's samples be: sample[i] /= (outSigArraySize / 2)
	const int16_t outSigArraySizeOverTwo = static_cast<int16_t>(outSigArraySize) / 2;
	for (size_t i = 0; i < componentArraySize; ++i)
	{
		copyOfDftRealComponent[i] /= outSigArraySizeOverTwo;
		copyOfDftComplexComponent[i] /= outSigArraySizeOverTwo;
	}

	// 0 the output array	
	const size_t outSizeMemSize = outSigArraySize * sizeof(int16_t);
	std::memset(outputSignal, 0, outSizeMemSize);

	// following are used for calculating the IDTF note that we're scaling from INT16 range to double range then convert back to store the calculated value
	const double iDftLength = static_cast<double>(outSigArraySize);
	double dblDftSampleIndex = 0.0, dblOutputSigSampleIndex = 0.0;

	// will need to put the int16 range values into float range values
	// double range (logical) = -1.0 to 1.0 (no capping)
	// INT16 range to -32768 t0 +32767
	constexpr double scaleToDoubleRange = 1.0 / 32767.0;
	constexpr double scaleToInt16Range = 32767.0 / 1.0;

	double realComponent, complexComponent;

	for (size_t dftSampleIndex = 0; dftSampleIndex < componentArraySize; ++dftSampleIndex, ++dblDftSampleIndex)
	{
		dblOutputSigSampleIndex = 0.0;
		for (size_t outputSignalSampleIndex = 0; outputSignalSampleIndex < outSigArraySize; ++outputSignalSampleIndex, ++dblOutputSigSampleIndex)
		{
			realComponent = static_cast<double>(copyOfDftRealComponent[dftSampleIndex]);
			complexComponent = static_cast<double>(copyOfDftComplexComponent[dftSampleIndex]);
			// scale value to double range
			realComponent *= scaleToDoubleRange;
			complexComponent *= scaleToDoubleRange;
			// perform the calculation
			realComponent *= cos(2.0 * M_PI * dblDftSampleIndex * dblOutputSigSampleIndex / iDftLength);
			complexComponent *= sin(2.0 * M_PI * dblDftSampleIndex * dblOutputSigSampleIndex / iDftLength);
			// scale value back into INT16 range
			realComponent *= scaleToInt16Range;
			complexComponent *= scaleToInt16Range;
			outputSignal[outputSignalSampleIndex] += static_cast<int16_t>(realComponent + complexComponent);
		}
	}

	// clean up don't want to leak memory
	free(copyOfDftRealComponent);
	free(copyOfDftComplexComponent);
}

void Signal::FourierTransforms::DiscreteFourierTransformFullF(const float* inputSignal, const size_t inputSignalLength, float* output)
{
	// 0 the output array
	const size_t outputArrayElements = OutputSamplesLength(inputSignalLength);
	const size_t outputBufferMemSz = sizeof(float) * outputArrayElements;
	std::memset(output, 0, outputBufferMemSz);

	// DiscreteFourierTransformF 
	float iAsFlt = 0.0f, jAsFlt = 0.0f;
	const float inputSignalLengthAsFlt = static_cast<float>(inputSignalLength);
	const float fltPi = static_cast<float>(M_PI);
	float currentFrequencyRealComponent = 0.0f, currentFrequencyComplexComponent = 0.0f;
	for (size_t i = 0; i < outputArrayElements; ++i, ++iAsFlt)
	{
		jAsFlt = 0.0f;
		currentFrequencyRealComponent = 0.0f, currentFrequencyComplexComponent = 0.0f;
		for (size_t j = 0; j < inputSignalLength; ++j, ++jAsFlt)
		{
			currentFrequencyRealComponent += inputSignal[j] * cosf(2.0f * fltPi * iAsFlt * jAsFlt / inputSignalLengthAsFlt);
			currentFrequencyComplexComponent -= inputSignal[j] * sinf(2.0f * fltPi * iAsFlt * jAsFlt / inputSignalLengthAsFlt);
		}
		output[i] = sqrtf((currentFrequencyRealComponent * currentFrequencyRealComponent)
			+ (currentFrequencyComplexComponent * currentFrequencyComplexComponent));
	}
}

void Signal::FourierTransforms::DiscreteFourierTransformFullTargetFrequenciesF(const float* inputSignal, const size_t inputSignalLength,
	const size_t numTargetFrequencies, const float* targetFrequencies, float* output)
{
	// 0 the output array
	const size_t outputBufferMemSz = sizeof(float) * numTargetFrequencies;
	std::memset(output, 0, outputBufferMemSz);

	float jAsFlt = 0.0f;
	const float inputSignalLengthAsFlt = static_cast<float>(inputSignalLength);
	const float fltPi = static_cast<float>(M_PI);
	float currentFrequencyRealComponent = 0.0f, currentFrequencyComplexComponent = 0.0f;
	for (size_t frequencyIndex = 0; frequencyIndex < numTargetFrequencies; ++frequencyIndex)
	{
		jAsFlt = 0.0f;
		currentFrequencyRealComponent = 0.0f, currentFrequencyComplexComponent = 0.0f;
		for (size_t j = 0; j < inputSignalLength; ++j, ++jAsFlt)
		{
			currentFrequencyRealComponent += inputSignal[j] * cosf(2.0f * fltPi * targetFrequencies[frequencyIndex] * jAsFlt / inputSignalLengthAsFlt);
			currentFrequencyComplexComponent -= inputSignal[j] * sinf(2.0f * fltPi * targetFrequencies[frequencyIndex] * jAsFlt / inputSignalLengthAsFlt);
		}
		output[frequencyIndex] = sqrtf((currentFrequencyRealComponent * currentFrequencyRealComponent)
			+ (currentFrequencyComplexComponent * currentFrequencyComplexComponent));
	}
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

void Signal::Filters::MovingAverageSubsquentPointsF(const float* inputSignal, const size_t inputSignalLength, float* output, size_t nPointsToAverage)
{
	const float nPointsToAvgAsFlt = static_cast<float>(nPointsToAverage);
	float average;

	for (size_t i = 0; i < inputSignalLength; ++i)
	{
		average = 0.0f;
		for (size_t j = i; j < inputSignalLength && j < nPointsToAverage + i; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsFlt;
		output[i] = average;
	}
}

void Signal::Filters::MovingAverageSubsquentPointsI16(const int16_t* inputSignal, const size_t inputSignalLength, int16_t* output, size_t nPointsToAverage)
{
	const int16_t nPointsToAvgAsI16 = static_cast<int16_t>(nPointsToAverage);
	int16_t average;

	for (size_t i = 0; i < inputSignalLength; ++i)
	{
		average = 0;
		for (size_t j = i; j < inputSignalLength && j < nPointsToAverage + i; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsI16;
		output[i] = average;
	}
}

void Signal::Filters::MovingAverageSymetricallyChosenPointsD(const double* inputSignal, const size_t inputSignalLength, double* output, size_t nPointsToAverage)
{
	const double nPointsToAvgAsDbl = static_cast<double>(nPointsToAverage);
	double average;
	const size_t halfNumPointsToAverage = nPointsToAverage / 2;

	for (int i = 0; i < inputSignalLength; ++i)
	{
		average = 0.0;
		const int posibleStartingIndex = i - static_cast<int>(halfNumPointsToAverage);
		const int sigAvgStartingPoint = posibleStartingIndex >= 0 ? posibleStartingIndex : 0; // start at i - halfNumPointsToAverage or the start of the array if we can't go back that far
		for (size_t j = sigAvgStartingPoint; j < inputSignalLength && j < sigAvgStartingPoint + nPointsToAverage; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsDbl;
		output[i] = average;
	}
}

void Signal::Filters::MovingAverageSymetricallyChosenPointsF(const float* inputSignal, const size_t inputSignalLength, float* output, size_t nPointsToAverage)
{
	const float nPointsToAvgAsFlt = static_cast<float>(nPointsToAverage);
	float average;

	const size_t halfNumPointsToAverage = nPointsToAverage / 2;

	for (int i = 0; i < inputSignalLength; ++i)
	{
		average = 0.0;
		const int posibleStartingIndex = i - static_cast<int>(halfNumPointsToAverage);
		const int sigAvgStartingPoint = posibleStartingIndex >= 0 ? posibleStartingIndex : 0; // start at i - halfNumPointsToAverage or the start of the array if we can't go back that far
		for (size_t j = sigAvgStartingPoint; j < inputSignalLength && j < sigAvgStartingPoint + nPointsToAverage; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsFlt;
		output[i] = average;
	}
}

void Signal::Filters::MovingAverageSymetricallyChosenPointsI16(const int16_t* inputSignal, const size_t inputSignalLength, int16_t* output, size_t nPointsToAverage)
{
	const int16_t nPointsToAvgAsFlt = static_cast<int16_t>(nPointsToAverage);
	int16_t average;
	const size_t halfNumPointsToAverage = nPointsToAverage / 2;

	for (int i = 0; i < inputSignalLength; ++i)
	{
		average = 0;
		const int posibleStartingIndex = i - static_cast<int>(halfNumPointsToAverage);
		const int sigAvgStartingPoint = posibleStartingIndex >= 0 ? posibleStartingIndex : 0; // start at i - halfNumPointsToAverage or the start of the array if we can't go back that far
		for (size_t j = sigAvgStartingPoint; j < inputSignalLength && j < sigAvgStartingPoint + nPointsToAverage; ++j)
		{
			average += inputSignal[j];
		}
		average /= nPointsToAvgAsFlt;
		output[i] = average;
	}
}

void Signal::Filters::RecursiveMovingAverageD(const double* intputSignal, const size_t inputSignalLength, double* outputSignal, const size_t nPointsToAverage)
{
	// 0 the output buffer
	const size_t outSigSZInMem = inputSignalLength * sizeof(double);
	memset(outputSignal, 0, outSigSZInMem);

	double accum = 0.0;
	for (size_t i = 0; i < nPointsToAverage - 1; ++i)
	{
		accum += intputSignal[i];
	}

	outputSignal[(nPointsToAverage - 1 / 2)] = accum / static_cast<double>(nPointsToAverage);
	const size_t halfNumPointsToAvg = nPointsToAverage / 2;
	for (size_t i = halfNumPointsToAvg; i < inputSignalLength - halfNumPointsToAvg - 1; ++i)
	{
		accum += intputSignal[i + ((nPointsToAverage - 1) / 2)] - intputSignal[i - halfNumPointsToAvg];
		outputSignal[i] = accum / static_cast<double>(nPointsToAverage);
	}
}

void Signal::Filters::RecursiveMovingAverageF(const float* intputSignal, const size_t inputSignalLength, float* outputSignal, const size_t nPointsToAverage)
{
	// 0 the output buffer
	const size_t outSigSZInMem = inputSignalLength * sizeof(float);
	memset(outputSignal, 0, outSigSZInMem);

	float accum = 0.0;
	for (size_t i = 0; i < nPointsToAverage - 1; ++i)
	{
		accum += intputSignal[i];
	}

	outputSignal[(nPointsToAverage - 1 / 2)] = accum / static_cast<float>(nPointsToAverage);
	const size_t halfNumPointsToAvg = nPointsToAverage / 2;
	for (size_t i = halfNumPointsToAvg; i < inputSignalLength - halfNumPointsToAvg - 1; ++i)
	{
		accum += intputSignal[i + ((nPointsToAverage - 1) / 2)] - intputSignal[i - halfNumPointsToAvg];
		outputSignal[i] = accum / static_cast<float>(nPointsToAverage);
	}
}

void Signal::Filters::RecursiveMovingAverageI16(const int16_t* intputSignal, const size_t inputSignalLength, int16_t* outputSignal, const size_t nPointsToAverage)
{
	// 0 the output buffer
	const size_t outSigSZInMem = inputSignalLength * sizeof(int16_t);
	memset(outputSignal, 0, outSigSZInMem);

	int16_t accum = 0;
	for (size_t i = 0; i < nPointsToAverage - 1; ++i)
	{
		accum += intputSignal[i];
	}

	outputSignal[(nPointsToAverage - 1 / 2)] = accum / static_cast<int16_t>(nPointsToAverage);
	const size_t halfNumPointsToAvg = nPointsToAverage / 2;
	for (size_t i = halfNumPointsToAvg; i < inputSignalLength - halfNumPointsToAvg - 1; ++i)
	{
		accum += intputSignal[i + ((nPointsToAverage - 1) / 2)] - intputSignal[i - halfNumPointsToAvg];
		outputSignal[i] = accum / static_cast<int16_t>(nPointsToAverage);
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

void Signal::Filters::Windowed::HammingF(float* outputWindow, const size_t inputSignalLength)
{
	// equation = w[i] = 0.54 - 0.46 * cos(2*PI*i / M)
	// i = sample index
	// M = number of samples

	const float fltInputSigLength = static_cast<float>(inputSignalLength);
	const float fltPi = static_cast<float>(M_PI); // done for type safety
	float iterAsFlt = 0.0f;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsFlt)
	{
		outputWindow[i] = 0.54f - 0.46f * cosf(2.0f * fltPi * iterAsFlt / fltInputSigLength);
	}
}

void Signal::Filters::Windowed::HammingI16(int16_t* outputWindow, const size_t inputSignalLength)
{
	// equation = w[i] = 0.54 - 0.46 * cos(2*PI*i / M)
	// i = sample index
	// M = number of samples

	const double dblInputSigLength = static_cast<double>(inputSignalLength);

	double iterAsDbl = 0.0;
	double valueAtI = 0.0;
	constexpr double scaleToInt16Range = 32767.0 / 1.0;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsDbl)
	{
		valueAtI = 0.54 - 0.46 * cos(2.0 * M_PI * iterAsDbl / dblInputSigLength);
		valueAtI *= scaleToInt16Range;
		outputWindow[i] = static_cast<int16_t>(valueAtI);
	}
}

void Signal::Filters::Windowed::BlackmanD(double* outputWindow, const size_t inputSignalLength)
{
	// equation = w[i] = 0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M)
	// i = sample index
	// M = number of samples

	const double dblInputSigLength = static_cast<double>(inputSignalLength);

	double iterAsDbl = 0.0;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsDbl)
	{
		outputWindow[i] = 0.42 - 0.5 * cos(2.0 * M_PI * iterAsDbl / dblInputSigLength) + 0.08 * cos(4.0 * M_PI * iterAsDbl / dblInputSigLength);
	}
}

void Signal::Filters::Windowed::BlackmanF(float* outputWindow, const size_t inputSignalLength)
{
	// equation = w[i] = 0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M)
	// i = sample index
	// M = number of samples

	const float fltInputSigLength = static_cast<float>(inputSignalLength);
	const float fltPi = static_cast<float>(M_PI); // done for type safety
	float iterAsFlt = 0.0f;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsFlt)
	{
		outputWindow[i] = 0.42f - 0.5f * cosf(2.0f * fltPi * iterAsFlt / fltInputSigLength) + 0.08f * cosf(4.0f * fltPi * iterAsFlt / fltInputSigLength);
	}
}

void Signal::Filters::Windowed::BlackmanI16(int16_t* outputWindow, const size_t inputSignalLength)
{
	// equation = w[i] = 0.42 - 0.5 * cos(2 * PI * i / M) + 0.08 * cos(4 * PI * i / M)
	// i = sample index
	// M = number of samples

	const double dblInputSigLength = static_cast<double>(inputSignalLength);
	double valueAtI = 0.0;
	constexpr double scaleToInt16Range = 32767.0 / 1.0;

	double iterAsDbl = 0.0;
	for (size_t i = 0; i < inputSignalLength; ++i, ++iterAsDbl)
	{
		valueAtI = 0.42 - 0.5 * cos(2.0 * M_PI * iterAsDbl / dblInputSigLength) + 0.08 * cos(4.0 * M_PI * iterAsDbl / dblInputSigLength);
		valueAtI *= scaleToInt16Range;
		outputWindow[i] = static_cast<int16_t>(valueAtI);
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

void Signal::Filters::Windowed::SyncLowPassF(float* filterOutput, const size_t filterOutputSize, float cutoffFrequency)
{
	const int filterOutputSizeI = static_cast<int>(filterOutputSize); // need signed values
	const int halfFilterLength = filterOutputSizeI / 2;

	float iFlt = 0.0f;
	const float filterOutputSizeFlt = static_cast<float>(filterOutputSize);
	const float fltPi = static_cast<float>(M_PI); // done for type safety
	for (int i = 0; i < filterOutputSizeI; ++i, ++iFlt)
	{
		if (i - halfFilterLength == 0) // half way point is a special value
		{
			filterOutput[i] = 2.0f * fltPi * cutoffFrequency;
		}
		else
		{
			filterOutput[i] = sinf(2.0f * fltPi * cutoffFrequency * (static_cast<float>(i - halfFilterLength)))
				/ (static_cast<float>(i - halfFilterLength));
			filterOutput[i] *= (0.54f - 0.46f * cosf(2.0f * fltPi * iFlt / filterOutputSizeFlt));
		}
	}
}

void Signal::Filters::Windowed::SyncLowPassI16(int16_t* filterOutput, const size_t filterOutputSize, double cutoffFrequency)
{
	const int filterOutputSizeI = static_cast<int>(filterOutputSize); // need signed values
	const int halfFilterLength = filterOutputSizeI / 2;

	double iDbl = 0.0;
	const double filterOutputSizeDbl = static_cast<double>(filterOutputSize);

	constexpr double scaleToInt16Range = 32767.0 / 1.0;
	double iThValue = 0.0;

	for (int i = 0; i < filterOutputSizeI; ++i, ++iDbl)
	{
		if (i - halfFilterLength == 0) // half way point is a special value
		{
			iThValue = 2.0 * M_PI * cutoffFrequency;
			iThValue *= scaleToInt16Range;
			filterOutput[i] = static_cast<int16_t>(iThValue);
		}
		else
		{
			iThValue = sin(2.0 * M_PI * cutoffFrequency * (static_cast<double>(i - halfFilterLength)))
				/ (static_cast<double>(i - halfFilterLength));
			iThValue *= (0.54 - 0.46 * cos(2.0 * M_PI * iDbl / filterOutputSizeDbl));
			iThValue *= scaleToInt16Range;
			filterOutput[i] = static_cast<int16_t>(iThValue);
		}
	}
}

void Signal::Normalisation::NormaliseF(float* signal, const size_t length)
{
	float max = 0.0f;
	for (size_t i = 0; i < length; ++i)
	{
		max = std::max(max, std::abs(signal[i]));
	}
	const float scaleBy = 1.0f / max;
	for (size_t i = 0; i < length; ++i)
	{
		signal[i] *= scaleBy;
	}
}

void Signal::Conversion::i16ToF(const int16_t* toConvert, const size_t length, float* output)
{
	// from (-32768 to 32767) to (-1.0f to 1.0f)
	static constexpr float SCALER = 1.0f / 32768.0f;
	for (size_t i = 0; i < length; ++i)
	{
		output[i] = static_cast<float>(toConvert[i]) * SCALER;
	}
}

void Signal::Conversion::fToI16(const float* toConvert, const size_t length, int16_t* output)
{
	// from (-1.0f to 1.0f) to (-32768 to 32767)
	static constexpr float SCALER = 32767.0f;
	for (size_t i = 0; i < length; ++i)
	{
		output[i] = static_cast<int16_t>(SCALER * toConvert[i]);
	}
}