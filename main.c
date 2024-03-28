#include <formatio.h>
#include <ansi_c.h>
#include <cvirte.h>		
#include <userint.h>
#include <utility.h>
#include <advanlys.h>
#include "main.h"
#include <analysis.h>
#include <toolbox.h>

#define SAMPLE_RATE		0
#define NPOINTS			1

//static int panelHandle;
static int panelFreq;
static int panelAcq;
int waveInfo[2];
double sampleRate = 0.0;
int numPoints = 0;
double* waveData = 0;
int* waveArray =0;
int fs=44100;
int startIndex =0;
int prevJpgIndex=1;
int nextJpgIndex=1;
double* py_envelope =0;

int main (int argc, char *argv[])
{   /*
	if (InitCVIRTE (0, argv, 0) == 0)
		return -1;	/* out of memory 
	if ((panelHandle = LoadPanel (0, "main.uir", PANEL)) < 0)
		return -1;
	DisplayPanel (panelHandle);
	RunUserInterface ();
	DiscardPanel (panelHandle);
	return 0;*/

	
	int error = 0;
	
	/* initialize and load resources */
	nullChk (InitCVIRTE (0, argv, 0));
	errChk (panelAcq = LoadPanel (0, "main.uir", PANEL));
	errChk (panelFreq = LoadPanel (0, "main.uir", PANEL_FREQ));
	
	/* display the panel and run the user interface */
	errChk (DisplayPanel (panelAcq));
	errChk (RunUserInterface ());

Error:
	/* clean up */
	if (panelAcq > 0)
		DiscardPanel (panelAcq);
	return 0;
}

int CVICALLBACK numberOfPassingThroughZero() {
	
	int number  =0;
	
	for(int i=1;i<numPoints;i++) {
		if(waveData[i]*waveData[i-1] <0)
			number++;
	}
	
	return number;
}

int CVICALLBACK loadWAVFile (int panel, int control, int event,
							 void *callbackData, int eventData1, int eventData2)
{
	
	switch (event)
	{
		case EVENT_COMMIT:
			
			FileToArray ("d:\\apd\\wafeInfo.txt", waveInfo, VAL_INTEGER, 2, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			sampleRate = waveInfo[SAMPLE_RATE];
			numPoints = waveInfo[NPOINTS];
			
			waveData = (double *) calloc(numPoints, sizeof(double));
			
			FileToArray ("d:\\apd\\waveData.txt", waveData, VAL_DOUBLE, numPoints, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			PlotY (panel, PANEL_GRAPH_RAW, waveData, numPoints, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			
			double maxValue = 0.0, minValue=0.0;
			int maxIndex =0, minIndex =0;
			double meanValue=0.0;
			double medianValue = 0.0;
			double dispersion =0.0;
			MaxMin1D(waveData, numPoints, &maxValue, &maxIndex, &minValue, &minIndex);
			Mean(waveData, numPoints, &meanValue);
			Median(waveData, numPoints, &medianValue);
			StdDev(waveData, numPoints,&meanValue, &dispersion);
			
			SetCtrlVal(panel, PANEL_MIN, minValue);
			SetCtrlVal(panel, PANEL_MAX, maxValue);
			SetCtrlVal(panel, PANEL_MEAN, meanValue);
			SetCtrlVal(panel, PANEL_MEDIAN, medianValue);
			SetCtrlVal(panel, PANEL_DISPERSION, dispersion);
			
			int numOfPassing = numberOfPassingThroughZero();
			SetCtrlVal(panel, PANEL_PASSING, numOfPassing);
			
			int intervals = 20;
			int* histogramArray = (int *) calloc (numPoints, sizeof(int));
			double* axisArray = (double *) calloc(intervals, sizeof(double));
			Histogram(waveData, numPoints,minValue-1, maxValue+1, histogramArray, axisArray, intervals);
			PlotXY(panel, PANEL_GRAPH_HISTOGRAM, axisArray, histogramArray,intervals,VAL_DOUBLE,VAL_INTEGER,VAL_VERTICAL_BAR,VAL_SOLID_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS,VAL_RED);
			
			double skewness, kurtosis;
			//skewness moment 3  + add control
			Moment(waveData, numPoints, 3, &skewness);
			//kurtosis moment 4
			Moment(waveData, numPoints, 3, &kurtosis);
			 
			SetCtrlVal(panel, PANEL_SKEWNESS, skewness);
			SetCtrlVal(panel, PANEL_KURTOSIS, kurtosis);
			
			int BitmapID;
			GetCtrlDisplayBitmap(panel, PANEL_GRAPH_RAW, 1, &BitmapID);
			SaveBitmapToJPEGFile(BitmapID, "D:\\apd\\pics\\graphFull.jpg", JPEG_PROGRESSIVE, 100);
			DiscardBitmap(BitmapID);
			/*
			FILE *myFile;
			int array_size;
    		myFile = fopen("d:\\apd\\envelope_size.txt", "r");
			fscanf(myFile, "%d", &array_size);
			
			py_envelope = (double *)calloc(array_size, sizeof(double));
			FileToArray ("d:\\apd\\amplitude_envelope.txt", py_envelope, VAL_DOUBLE, array_size, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			
			PlotY (panel, PANEL_GRAPH_RAW, py_envelope, array_size, VAL_DOUBLE, VAL_SCATTER, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_GREEN);
			*/
			
			//anvelopa
			int count;
			double* peakLocations = (double *)malloc(numPoints * sizeof(double));
			double* peakAmplitudes = (double *)malloc(numPoints * sizeof(double));
			double* peakSecondDerivatives = (double *)malloc(numPoints * sizeof(double));
			
			PeakDetector(waveData, numPoints, 100,3,0,1,1,&count, &peakLocations, &peakAmplitudes, &peakSecondDerivatives);
		
			double* vmax = (double *)malloc(numPoints * sizeof(double));
			int* indmax = (int *)malloc(numPoints * sizeof(int));
			int k=0;
			for(int i=0;i<count;i+=50) { 
				
				double max =peakAmplitudes[i];  
			    int maxind = peakLocations[i];
			    for(int j=i;j<i+50;j++) {
			        if(peakAmplitudes[j] > max){ 
						max = peakAmplitudes[j];
						maxind = peakLocations[j];
					}
				}
				vmax[k] = max;
				indmax[k] = maxind;
				k++;
			
			}
			PlotXY(panel, PANEL_GRAPH_RAW, indmax, vmax, k-1, VAL_INTEGER, VAL_DOUBLE, VAL_FAT_LINE, VAL_SOLID_CIRCLE, VAL_DOT,  VAL_CONNECTED_POINTS, VAL_GREEN);
			
			
			break;
	}
	return 0;
}

int CVICALLBACK goPrevInterval (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	int* inter = (int *) calloc(fs, sizeof(int));
	switch (event)
	{
		case EVENT_COMMIT:
			
			if(startIndex >= fs && startIndex <= fs*6-1) {
				startIndex -= fs;
				for(int i=0;i<fs;i++)
					inter[i] = i+startIndex;
				DeleteGraphPlot(panel, PANEL_GRAPH_RAW, -1, VAL_IMMEDIATE_DRAW);
				
				PlotXY(panel, PANEL_GRAPH_RAW,inter ,waveData + startIndex, fs, VAL_INTEGER, VAL_DOUBLE,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
				
				SetCtrlVal(panel, PANEL_NUMERIC_START, startIndex/fs);
				SetCtrlVal(panel, PANEL_NUMERIC_STOP, (startIndex+fs)/fs);
				
				//save picture of second
				int BitmapID;
				char path[28] = "D:\\apd\\pics\\graphPrev";
				char tmp = prevJpgIndex + '0';
				char* jpg = (char *)(strncat(path, &tmp, 1));
				jpg = strcat(jpg, ".jpg");
				
				GetCtrlDisplayBitmap(panel, PANEL_GRAPH_RAW, 1, &BitmapID);
				SaveBitmapToJPEGFile(BitmapID, jpg, JPEG_PROGRESSIVE, 100);
				DiscardBitmap(BitmapID);
				
				prevJpgIndex++;
			}
			break;
	}
	return 0;
}

int CVICALLBACK goNextInterval (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	int* inter = (int *) calloc(fs, sizeof(int));
	
	switch (event)
	{
		case EVENT_COMMIT:
			
			if(startIndex >= 0 &&  startIndex <= fs*5-1) {
				
				startIndex += fs;
				for(int i=0;i<fs;i++)
					inter[i] = i+startIndex;
				DeleteGraphPlot(panel, PANEL_GRAPH_RAW, -1, VAL_IMMEDIATE_DRAW);
				
				PlotXY(panel, PANEL_GRAPH_RAW,inter ,waveData + startIndex, fs, VAL_INTEGER, VAL_DOUBLE,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
				SetCtrlVal(panel, PANEL_NUMERIC_START, startIndex/fs);
				SetCtrlVal(panel, PANEL_NUMERIC_STOP, (startIndex+fs)/fs);
				
				//save picture of second
				int BitmapID;
				char path[28] = "D:\\apd\\pics\\graphNext";
				char tmp = nextJpgIndex + '0';
				char* jpg = (char *)(strncat(path, &tmp, 1));
				jpg = strcat(jpg, ".jpg");
				
				GetCtrlDisplayBitmap(panel, PANEL_GRAPH_RAW, 1, &BitmapID);
				SaveBitmapToJPEGFile(BitmapID, jpg, JPEG_PROGRESSIVE, 100);
				DiscardBitmap(BitmapID);
				nextJpgIndex++;
			}
			break;
	}
	return 0;
}

int CVICALLBACK applyFilter (int panel, int control, int event,
							 void *callbackData, int eventData1, int eventData2)
{
	double alpha=0;
	int option, start, stop, window,dim;
	int* inter=0;
	
	switch (event)
	{
		case EVENT_COMMIT:
			
			GetCtrlVal(panel, PANEL_RING_FILTER, &option);
			GetCtrlVal(panel, PANEL_NUMERIC_START, &start);
			GetCtrlVal(panel, PANEL_NUMERIC_STOP, &stop);
			GetCtrlVal(panel, PANEL_WINDOW, &window);
			DeleteGraphPlot(panel, PANEL_GRAPH_FILTERED, -1, VAL_IMMEDIATE_DRAW);
		
			if(option ==0) //mediere pe 32
			{
				if(start==stop) 
					dim = numPoints;
				else dim = fs;
				
				double* medSum = (double *) calloc(dim, sizeof(double));
				inter = (int *) calloc(dim, sizeof(int));
                double sum=0;
				int z=0;
				
				for(int j=0;j<window;j++) 
						sum+= waveData[j+startIndex];	
				medSum[z++] = (double)sum/window;
				
				for(int i=1;i<dim-window;i++)
					
				{   
					sum = sum + waveData[window+startIndex+i-1] - waveData[startIndex+i-1] ;
					medSum[z] = (double)sum/window;
				    z++;
					
				}
				for(int i=dim-window;i<dim;i++)	
				{
					medSum[z] = medSum[z-1];
					z++;
				}
				for(int i=0;i<dim;i++) {
						inter[i] = i+startIndex;
					}
				PlotXY(panel, PANEL_GRAPH_FILTERED, inter, medSum,dim, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
				
			}
			else if(option==1) //de ordin 1
			{
				GetCtrlVal(panel, PANEL_ALPHA, &alpha);
				
				if(start==stop) 
					dim = numPoints;
				else dim = fs;
				
				double* filt = (double *) calloc(dim, sizeof(double));
				inter = (int *) calloc(dim, sizeof(int));
				filt[0]=0;
				
				for(int i=1;i<dim;++i) {
						filt[i] = (1-alpha)*filt[i-1]+alpha*waveData[i+startIndex];
				}
				
				for(int i=0;i<dim;i++) 
					inter[i] = i+startIndex;
				
				PlotXY(panel, PANEL_GRAPH_FILTERED, inter, filt,dim, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
				
				break;
			}
			else if(option==2) {
				//derivata
				
				if(start==stop) 
					dim = numPoints;
				else dim = fs;
				
				double* deriv = (double *) calloc(dim, sizeof(double));
				inter = (int *) calloc(dim, sizeof(int));
				
				deriv[0]=0;
				for(int i=1;i<dim;i++)
					deriv[i] = waveData[i]-waveData[i-1];
				
				for(int i=0;i<dim;i++) 
					inter[i] = i+startIndex;
				
				PlotXY(panel, PANEL_GRAPH_FILTERED, inter, deriv, dim, VAL_INTEGER, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
				break;
			}
			break;
	}
	return 0;
}



int CVICALLBACK onMainPanel (int panel, int event, void *callbackData,
							 int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:

			break;
		case EVENT_LOST_FOCUS:

			break;
		case EVENT_CLOSE:
			QuitUserInterface(0);
			break;
	}
	return 0;
}

int CVICALLBACK onSwitch (int panel, int control, int event,
						  void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			if(panel == panelAcq)
			{
				SetCtrlVal(panelFreq, PANEL_FREQ_BINARYSWITCH, 1);
				DisplayPanel(panelFreq);
				HidePanel(panel);
			}
			else
			{
				SetCtrlVal(panelAcq, PANEL_FREQ_BINARYSWITCH, 0);
				DisplayPanel(panelAcq);
				HidePanel(panel);
			}
			
			break;
	}
	return 0;
}

int CVICALLBACK OnTimer (int panel, int control, int event,
						 void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_TIMER_TICK:

			break;
	}
	return 0;
}

double* CVICALLBACK downsampleSpectrum(double* spectrum, int n) {
	
	double* downsampled = (double *)malloc(n/2*sizeof(double));
	for(int i=0;i<n;i+=2)
	{
		downsampled[i/2] = (spectrum[i] + spectrum[i+1])/2;
	}
	return downsampled;
}

int CVICALLBACK onLoadSpectrum (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	
	WindowConst winConst;
	int n,second;
	double df,frequencyPeak, powerPeak;
	
	char unit[32] = "V";
	
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelFreq,PANEL_FREQ_RING_N, &n);
			GetCtrlVal(panelFreq, PANEL_FREQ_NUMERIC_SECOND, &second);
			double* spectrumArray = (double *)malloc(n*sizeof(double));
			double* autoSpectrum = (double *)malloc(n/2*sizeof(double));
			double* convertedSpectrum = (double *)malloc(n/2*sizeof(double));
			double* downsampled;
			int tn = n;
			if(second==0) {
				
				Copy1D(waveData, n, spectrumArray);
				downsampled = downsampleSpectrum(spectrumArray, n);
				//n=n/2;
				Copy1D(downsampled, n/2, spectrumArray);
				ScaledWindowEx (spectrumArray, n, RECTANGLE, 0 ,&winConst);
				AutoPowerSpectrum(spectrumArray, n, 1.0/(sampleRate/2), autoSpectrum, &df);
				PowerFrequencyEstimate (autoSpectrum, n/2, -1, winConst, df, 7, &frequencyPeak, &powerPeak);
				SpectrumUnitConversion (autoSpectrum, n/2, 0, 0, 0, df, winConst, convertedSpectrum, unit);
				DeleteGraphPlot(panelFreq, PANEL_FREQ_GRAPH_SPECTRUM, -1, VAL_IMMEDIATE_DRAW);
				
				PlotWaveform(panelFreq, PANEL_FREQ_GRAPH_SPECTRUM, convertedSpectrum, n/2, VAL_DOUBLE, 1,0,0,df,VAL_THIN_LINE, VAL_SOLID_SQUARE,
							 VAL_SOLID, VAL_CONNECTED_POINTS,VAL_RED);
				
			}		 
			else {
				
				Copy1D(waveData+fs*second, n, spectrumArray);
				downsampled = downsampleSpectrum(spectrumArray, n);
				//n=n/2;
				Copy1D(downsampled, n/2, spectrumArray);
				ScaledWindowEx (spectrumArray, n, RECTANGLE, 0 ,&winConst);
				AutoPowerSpectrum(spectrumArray, n, 1.0/(sampleRate/2), autoSpectrum, &df);
				PowerFrequencyEstimate (autoSpectrum, n/2, -1, winConst, df, 7, &frequencyPeak, &powerPeak);
				SpectrumUnitConversion (autoSpectrum, n/2, 0, 0, 0, df, winConst, convertedSpectrum, unit);
				DeleteGraphPlot(panelFreq, PANEL_FREQ_GRAPH_SPECTRUM, -1, VAL_IMMEDIATE_DRAW);
				PlotWaveform(panelFreq, PANEL_FREQ_GRAPH_SPECTRUM, convertedSpectrum, n/2, VAL_DOUBLE, 1,0,0,df,VAL_THIN_LINE, VAL_SOLID_SQUARE,
							 VAL_SOLID, VAL_CONNECTED_POINTS,VAL_RED);
			}
			
			if(second==0) {
				DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_RAW,-1,VAL_IMMEDIATE_DRAW); 
				PlotY (panel, PANEL_FREQ_GRAPH_RAW, waveData, tn, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			}
			else {
				DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_RAW,-1,VAL_IMMEDIATE_DRAW); 
				PlotY (panel, PANEL_FREQ_GRAPH_RAW, waveData+fs*second, tn, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			}
			
			int BitmapID, month,day,year,h,m,s;
			char path[91] = "D:\\apd\\pics\\signalBeforeFiltering";
			char* custom=(char *)malloc(sizeof(char)*46);
			GetSystemDate(&month, &day, &year);
			GetSystemTime(&h, &m,&s);
			
			Fmt(custom, "%d_%d_%d__%d_%d_%d__second_%d__npoints__%d", day,month,year,h,m,s,second, n);
			char* jpg = (char *)(strcat(path, custom));
			jpg = strcat(jpg, ".jpg");
			GetCtrlDisplayBitmap(panel, PANEL_FREQ_GRAPH_RAW, 1, &BitmapID);
			SaveBitmapToJPEGFile(BitmapID, jpg, JPEG_PROGRESSIVE, 100);
			
			
			strcpy(path,"D:\\apd\\pics\\spectrumBeforeFiltering");
			Fmt(custom, "%d_%d_%d__%d_%d_%d__second_%d__npoints__%d", day,month,year,h,m,s,second, n);
			jpg = (char *)(strcat(path, custom));
			jpg = strcat(jpg, ".jpg");
			GetCtrlDisplayBitmap(panel, PANEL_FREQ_GRAPH_SPECTRUM, 1, &BitmapID);
			SaveBitmapToJPEGFile(BitmapID, jpg, JPEG_PROGRESSIVE, 100);
			
			DiscardBitmap(BitmapID);
			
			break;
	}
	return 0;
}
///SS LA SPECTRU
int CVICALLBACK onFilterButton (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	int filterType,second,vs;
	int windowType,n;
	WindowConst winConst;
	switch (event)
	{
		case EVENT_COMMIT:
			GetCtrlVal(panelFreq, PANEL_FREQ_RING_FILTER, &filterType);
			GetCtrlVal(panelFreq, PANEL_FREQ_RING_WINDOW, &windowType);
			GetCtrlVal(panelFreq, PANEL_FREQ_RING_N, &n);
			GetCtrlVal(panelFreq, PANEL_FREQ_NUMERIC_SECOND, &second);
			double* filteredData;
			double* tempData = (double *)malloc(sizeof(double)*n);
			
			if(second ==0) {
			
				Copy1D(waveData, n, tempData);
			}
			else {
				
				Copy1D(waveData+fs*second, n, tempData);
			}
			
			DeleteGraphPlot(panel, PANEL_FREQ_GRAPH_FILTERED,-1,VAL_IMMEDIATE_DRAW); 
			
			if(filterType==0) {
				//kaiser
				//Consideram sampling frequency = fs, cutoff frequency = 900, filter length = 55, beta = 4.5
				double fcut = 900, beta=4.5;
				int  fl =55;
				double* coefArray = (double *)malloc(sizeof(double)*fl);
				
				Ksr_HPF(fs/2,fcut,fl, coefArray, beta);
				vs=n+fl-1;
				filteredData =  (double *)malloc(sizeof(double)*vs);
				Convolve(coefArray, fl, tempData, n, filteredData);
				 
				
				
			} else if(filterType==1) {
				//butterworth
				double fcut=900, fpass=700, fsampling=fs; // fc<=fs/2 
				vs=n;
				filteredData =  (double *)malloc(sizeof(double)*vs);
				Bw_HPF(tempData, n,fsampling/2 , fcut, 6, filteredData); 
				
			
			}
			double* tempFilteredData = (double *)malloc(sizeof(double)*vs);
			
			if(windowType==0) {
				//Hamming
				ScaledWindowEx(filteredData, vs, HAMMING, 0 ,&winConst);
				
			}
			else  if(windowType==1) {
				//Welch
				ScaledWindowEx(filteredData, vs, WELCH, 0 ,&winConst);
			}
			
			PlotY (panel, PANEL_FREQ_GRAPH_FILTERED, filteredData, vs, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			
			double* spectrumArray = (double *)malloc(n*sizeof(double));
			double* autoSpectrum = (double *)malloc(n/2*sizeof(double));
			double* convertedSpectrum = (double *)malloc(n/2*sizeof(double));
			double df,frequencyPeak, powerPeak;
			double* downsampled;
	        char unit[32] = "V";
			
			if(second==0) {
		
				Copy1D(filteredData, n, spectrumArray);
				downsampled = downsampleSpectrum(spectrumArray, n);
				//n=n/2;
				
				Copy1D(downsampled, n/2, spectrumArray);
				ScaledWindowEx (spectrumArray, n, RECTANGLE, 0 ,&winConst);
				AutoPowerSpectrum(spectrumArray, n, 1.0/(sampleRate/2), autoSpectrum, &df);
				PowerFrequencyEstimate (autoSpectrum, n/2, -1, winConst, df, 7, &frequencyPeak, &powerPeak);
				SpectrumUnitConversion (autoSpectrum, n/2, 0, 0, 0, df, winConst, convertedSpectrum, unit);
				DeleteGraphPlot(panelFreq, PANEL_FREQ_GRAPH_FILT_SPECTRUM, -1, VAL_IMMEDIATE_DRAW);
				PlotWaveform(panelFreq, PANEL_FREQ_GRAPH_FILT_SPECTRUM, convertedSpectrum, n/2, VAL_DOUBLE, 1,0,0,df,VAL_THIN_LINE, VAL_SOLID_SQUARE,
							 VAL_SOLID, VAL_CONNECTED_POINTS,VAL_RED);
			}		 
			else {
				
				Copy1D(filteredData, n, spectrumArray);
				downsampled = downsampleSpectrum(spectrumArray, n);
				//n=n/2;
				
				Copy1D(downsampled, n/2, spectrumArray);
				ScaledWindowEx (spectrumArray, n, RECTANGLE, 0 ,&winConst);
				AutoPowerSpectrum(spectrumArray, n, 1.0/(sampleRate/2), autoSpectrum, &df);
				PowerFrequencyEstimate (autoSpectrum, n/2, -1, winConst, df, 7, &frequencyPeak, &powerPeak);
				SpectrumUnitConversion (autoSpectrum, n/2, 0, 0, 0, df, winConst, convertedSpectrum, unit);
				DeleteGraphPlot(panelFreq, PANEL_FREQ_GRAPH_FILT_SPECTRUM, -1, VAL_IMMEDIATE_DRAW);
				PlotWaveform(panelFreq, PANEL_FREQ_GRAPH_FILT_SPECTRUM, convertedSpectrum, n/2, VAL_DOUBLE, 1,0,0,df,VAL_THIN_LINE, VAL_SOLID_SQUARE,
							 VAL_SOLID, VAL_CONNECTED_POINTS,VAL_RED);
			}
			
			//ArrayToFile("original.txt", waveData, VAL_DOUBLE, vs, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_SEP_BY_TAB,0, VAL_ASCII, VAL_TRUNCATE);
			//ArrayToFile("filtered.txt", filteredData, VAL_DOUBLE,vs, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_SEP_BY_TAB,0, VAL_ASCII, VAL_TRUNCATE); 
			
			int BitmapID, month,day,year,h,m,s;
			char path[91] = "D:\\apd\\pics\\signalAfterFiltering";
			char* custom=(char *)malloc(sizeof(char)*46);
			GetSystemDate(&month, &day, &year);
			GetSystemTime(&h, &m,&s);
			
			Fmt(custom, "%d_%d_%d__%d_%d_%d__second_%d__npoints__%d", day,month,year,h,m,s,second, n);
			char* jpg = (char *)(strcat(path, custom));
			jpg = strcat(jpg, ".jpg");
			GetCtrlDisplayBitmap(panel, PANEL_FREQ_GRAPH_FILTERED, 1, &BitmapID);
			SaveBitmapToJPEGFile(BitmapID, jpg, JPEG_PROGRESSIVE, 100);
			
			
			strcpy(path,"D:\\apd\\pics\\spectrumAfterFiltering");
			Fmt(custom, "%d_%d_%d__%d_%d_%d__second_%d__npoints__%d", day,month,year,h,m,s,second, n);
			jpg = (char *)(strcat(path, custom));
			jpg = strcat(jpg, ".jpg");
			GetCtrlDisplayBitmap(panel, PANEL_FREQ_GRAPH_FILT_SPECTRUM, 1, &BitmapID);
			SaveBitmapToJPEGFile(BitmapID, jpg, JPEG_PROGRESSIVE, 100);
			
			DiscardBitmap(BitmapID);
			break;
	}
	return 0;
}


int CVICALLBACK onFreqPanel (int panel, int event, void *callbackData,
							 int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:

			break;
		case EVENT_LOST_FOCUS:

			break;
		case EVENT_CLOSE:
			QuitUserInterface(0);
			break;
	}
	return 0;
}

int CVICALLBACK onSecondSwitch (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			if(panel == panelAcq)
			{
				SetCtrlVal(panelFreq, PANEL_BINARYSWITCH, 1);
				DisplayPanel(panelFreq);
				HidePanel(panel);
			}
			else
			{
				SetCtrlVal(panelAcq, PANEL_BINARYSWITCH, 0);
				DisplayPanel(panelAcq);
				HidePanel(panel);
			}
			
			break;
	}
	return 0;
}
