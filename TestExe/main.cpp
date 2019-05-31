#include <iostream>
#include <fstream>

#include <Signal.h>

namespace testData
{
	// below array(s) are taken from the course's example code
	/* ----------------------------------------------------------------------
	** Test input signal contains 1000Hz + 15000 Hz
	** ------------------------------------------------------------------- */

	double InputSignal_f32_1kHz_15kHz[320] =
	{
	+0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	-0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	-0.8660254038f, -0.4619397663f, -1.3194792169f, -1.1827865776f, -0.5000000000f, -1.1827865776f, -1.3194792169f, -0.4619397663f,
	-0.8660254038f, -1.2552931065f, -0.3535533906f, -0.4174197128f, -1.0000000000f, -0.1913417162f, +0.0947343455f, -0.5924659585f,
	+0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	+0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	-0.8660254038f, -0.4619397663f, -1.3194792169f, -1.1827865776f, -0.5000000000f, -1.1827865776f, -1.3194792169f, -0.4619397663f,
	-0.8660254038f, -1.2552931065f, -0.3535533906f, -0.4174197128f, -1.0000000000f, -0.1913417162f, +0.0947343455f, -0.5924659585f,
	+0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	+0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	-0.8660254038f, -0.4619397663f, -1.3194792169f, -1.1827865776f, -0.5000000000f, -1.1827865776f, -1.3194792169f, -0.4619397663f,
	-0.8660254038f, -1.2552931065f, -0.3535533906f, -0.4174197128f, -1.0000000000f, -0.1913417162f, +0.0947343455f, -0.5924659585f,
	-0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	-0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	-0.8660254038f, -0.4619397663f, -1.3194792169f, -1.1827865776f, -0.5000000000f, -1.1827865776f, -1.3194792169f, -0.4619397663f,
	-0.8660254038f, -1.2552931065f, -0.3535533906f, -0.4174197128f, -1.0000000000f, -0.1913417162f, +0.0947343455f, -0.5924659585f,
	+0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	+0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	-0.8660254038f, -0.4619397663f, -1.3194792169f, -1.1827865776f, -0.5000000000f, -1.1827865776f, -1.3194792169f, -0.4619397663f,
	-0.8660254038f, -1.2552931065f, -0.3535533906f, -0.4174197128f, -1.0000000000f, -0.1913417162f, +0.0947343455f, -0.5924659585f,
	-0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	+0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	-0.8660254038f, -0.4619397663f, -1.3194792169f, -1.1827865776f, -0.5000000000f, -1.1827865776f, -1.3194792169f, -0.4619397663f,
	-0.8660254038f, -1.2552931065f, -0.3535533906f, -0.4174197128f, -1.0000000000f, -0.1913417162f, +0.0947343455f, -0.5924659585f,
	-0.0000000000f, +0.5924659585f, -0.0947343455f, +0.1913417162f, +1.0000000000f, +0.4174197128f, +0.3535533906f, +1.2552931065f,
	+0.8660254038f, +0.4619397663f, +1.3194792169f, +1.1827865776f, +0.5000000000f, +1.1827865776f, +1.3194792169f, +0.4619397663f,
	+0.8660254038f, +1.2552931065f, +0.3535533906f, +0.4174197128f, +1.0000000000f, +0.1913417162f, -0.0947343455f, +0.5924659585f,
	+0.0000000000f, -0.5924659585f, +0.0947343455f, -0.1913417162f, -1.0000000000f, -0.4174197128f, -0.3535533906f, -1.2552931065f,
	};

	//Low-pass filter with 6000Hz cutoff frequency
	double  Impulse_response[29] = {
	 -0.0018225230f, -0.0015879294f, +0.0000000000f, +0.0036977508f, +0.0080754303f, +0.0085302217f, -0.0000000000f, -0.0173976984f,
	 -0.0341458607f, -0.0333591565f, +0.0000000000f, +0.0676308395f, +0.1522061835f, +0.2229246956f, +0.2504960933f, +0.2229246956f,
	 +0.1522061835f, +0.0676308395f, +0.0000000000f, -0.0333591565f, -0.0341458607f, -0.0173976984f, -0.0000000000f, +0.0085302217f,
	 +0.0080754303f, +0.0036977508f, +0.0000000000f, -0.0015879294f, -0.0018225230f
	};
}

void DumpWaveformToTextFileD(const char* filePath, const double* waveformSamlpes, const size_t nWaveFormSamples)
{
	std::ofstream output(filePath);
	for (size_t i = 0; i < nWaveFormSamples; ++i)
	{
		output << waveformSamlpes[i] << std::endl;
	}
	output.close();
}

int main()
{
	const size_t inputSignalArrayElementCount = sizeof(testData::InputSignal_f32_1kHz_15kHz) / sizeof(testData::InputSignal_f32_1kHz_15kHz[0]);
	const double waveFormMean = Signal::Statistics::MeanD(testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount);
	const double waveFormVarience = Signal::Statistics::VarienceD(testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount, waveFormMean);
	const double waveFormStdDeviation = Signal::Statistics::StandardDeviationD(waveFormVarience);
	std::cout << "Waveform mean calculated as: " << waveFormMean << std::endl;
	std::cout << "Waveform varience calculated as: " << waveFormVarience << std::endl;
	std::cout << "Waveform standard deviation calculated as: " << waveFormStdDeviation << std::endl;

	// convolution parts of the library
	const size_t impulseResArrayMemSize = sizeof(testData::Impulse_response);
	const size_t impulseResArrayElementCount = impulseResArrayMemSize / sizeof(testData::Impulse_response[0]);
	const size_t outSignalSize = inputSignalArrayElementCount + impulseResArrayElementCount;
	double outSignal[outSignalSize];
	Signal::Convolution::ConvolutionD(testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount, testData::Impulse_response, impulseResArrayElementCount, outSignal);

	double runningSumArray[inputSignalArrayElementCount];
	Signal::Convolution::RunningSumD(testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount, runningSumArray);

	double convDifferenceArray[inputSignalArrayElementCount];
	Signal::Convolution::DifferenceD(testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount, convDifferenceArray);

	// DFT part of the library
	const size_t dftOutputArrayElementCount = (inputSignalArrayElementCount / 2) + 1;
	double waveformRealComponent[dftOutputArrayElementCount];
	double waveformComplexComponent[dftOutputArrayElementCount];
	double dftMag[dftOutputArrayElementCount];
	Signal::FourierTransforms::DiscreteFourierTransformD(testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount, waveformRealComponent, waveformComplexComponent);
	Signal::FourierTransforms::DiscreteFourierTransformMagnitudeD(dftMag, waveformRealComponent, waveformComplexComponent, dftOutputArrayElementCount);
	double idftArray[inputSignalArrayElementCount];
	Signal::FourierTransforms::InverseDiscreteFourierTransformD(idftArray, waveformRealComponent, waveformComplexComponent, dftOutputArrayElementCount);

	// write signals to file (human readable)
	DumpWaveformToTextFileD("ConvolutedSignal.Signal", outSignal, outSignalSize);
	DumpWaveformToTextFileD("InputSignal.Signal", testData::InputSignal_f32_1kHz_15kHz, inputSignalArrayElementCount);
	DumpWaveformToTextFileD("ImpulseResponce.Signal", testData::Impulse_response, impulseResArrayElementCount);
	DumpWaveformToTextFileD("RunningSum.Signal", runningSumArray, inputSignalArrayElementCount);
	DumpWaveformToTextFileD("ConvolutionDifference.Signal", convDifferenceArray, inputSignalArrayElementCount);
	DumpWaveformToTextFileD("DFT_RealComponent.Signal", waveformRealComponent, dftOutputArrayElementCount);
	DumpWaveformToTextFileD("DFT_ComplexComponent.Signal", waveformRealComponent, dftOutputArrayElementCount);
	DumpWaveformToTextFileD("DFT_Magnitude.Signal", dftMag, dftOutputArrayElementCount);
	DumpWaveformToTextFileD("IDFT.Signal", idftArray, inputSignalArrayElementCount);

	return 0;
}