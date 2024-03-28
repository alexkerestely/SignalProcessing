/**************************************************************************/
/* LabWindows/CVI User Interface Resource (UIR) Include File              */
/*                                                                        */
/* WARNING: Do not add to, delete from, or otherwise modify the contents  */
/*          of this include file.                                         */
/**************************************************************************/

#include <userint.h>

#ifdef __cplusplus
    extern "C" {
#endif

     /* Panels and Controls: */

#define  PANEL                            1       /* callback function: onMainPanel */
#define  PANEL_GRAPH_FILTERED             2       /* control type: graph, callback function: (none) */
#define  PANEL_GRAPH_RAW                  3       /* control type: graph, callback function: (none) */
#define  PANEL_BUTTON_APPLY               4       /* control type: command, callback function: applyFilter */
#define  PANEL_BUTTON_NEXT                5       /* control type: command, callback function: goNextInterval */
#define  PANEL_BUTTON_PREV                6       /* control type: command, callback function: goPrevInterval */
#define  PANEL_RING_FILTER                7       /* control type: ring, callback function: (none) */
#define  PANEL_BUTTON_LOAD                8       /* control type: command, callback function: loadWAVFile */
#define  PANEL_PASSING                    9       /* control type: numeric, callback function: (none) */
#define  PANEL_DISPERSION                 10      /* control type: numeric, callback function: (none) */
#define  PANEL_MEDIAN                     11      /* control type: numeric, callback function: (none) */
#define  PANEL_MEAN                       12      /* control type: numeric, callback function: (none) */
#define  PANEL_MAX                        13      /* control type: numeric, callback function: (none) */
#define  PANEL_KURTOSIS                   14      /* control type: numeric, callback function: (none) */
#define  PANEL_SKEWNESS                   15      /* control type: numeric, callback function: (none) */
#define  PANEL_MIN                        16      /* control type: numeric, callback function: (none) */
#define  PANEL_GRAPH_HISTOGRAM            17      /* control type: graph, callback function: (none) */
#define  PANEL_TEXTMSG                    18      /* control type: textMsg, callback function: (none) */
#define  PANEL_NUMERIC_STOP               19      /* control type: numeric, callback function: (none) */
#define  PANEL_WINDOW                     20      /* control type: numeric, callback function: (none) */
#define  PANEL_NUMERIC_START              21      /* control type: numeric, callback function: (none) */
#define  PANEL_ALPHA                      22      /* control type: numeric, callback function: (none) */
#define  PANEL_BINARYSWITCH               23      /* control type: binary, callback function: onSwitch */

#define  PANEL_FREQ                       2       /* callback function: onFreqPanel */
#define  PANEL_FREQ_BINARYSWITCH          2       /* control type: binary, callback function: onSecondSwitch */
#define  PANEL_FREQ_GRAPH_SPECTRUM        3       /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_LOAD_SPECTRUM_BUTTON  4       /* control type: command, callback function: onLoadSpectrum */
#define  PANEL_FREQ_RING_N                5       /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_NUMERIC_SECOND        6       /* control type: numeric, callback function: (none) */
#define  PANEL_FREQ_RING_WINDOW           7       /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_RING_FILTER           8       /* control type: ring, callback function: (none) */
#define  PANEL_FREQ_FILTER_BUTTON         9       /* control type: command, callback function: onFilterButton */
#define  PANEL_FREQ_GRAPH_RAW             10      /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH_FILTERED        11      /* control type: graph, callback function: (none) */
#define  PANEL_FREQ_GRAPH_FILT_SPECTRUM   12      /* control type: graph, callback function: (none) */


     /* Control Arrays: */

          /* (no control arrays in the resource file) */


     /* Menu Bars, Menus, and Menu Items: */

          /* (no menu bars in the resource file) */


     /* Callback Prototypes: */

int  CVICALLBACK applyFilter(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK goNextInterval(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK goPrevInterval(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK loadWAVFile(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK onFilterButton(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK onFreqPanel(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK onLoadSpectrum(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK onMainPanel(int panel, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK onSecondSwitch(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);
int  CVICALLBACK onSwitch(int panel, int control, int event, void *callbackData, int eventData1, int eventData2);


#ifdef __cplusplus
    }
#endif