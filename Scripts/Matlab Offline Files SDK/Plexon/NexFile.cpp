/* 

  1. NEX FILE STRUCTURE
  ------------------------

  Nex data file has the following structure:

	- file header (structure NexFileHeader)
	- variable headers (structures NexVarHeader)
	- variable data (stored as binary arrays, depends on data type)

  Each variable header contains the size of the array that stores
  the variable data as well as the location of this array in the file.

  Here is the minimal code that reads Nex data file:

  -----------------------------------------------------
  void ReadNex(){
	NexFileHeader fh;
	// assume that file has less than 1000 variables
	NexVarHeader  vh[1000];
	int i;
	FILE* fp;

	fp = fopen("text.nex", "rb");
	// 1. read file header
	fread(&fh, sizeof(NexFileHeader), 1, fp);

   // 2. read the variable headers
	for(i=0; i<fh.NumVars; i++){
		fread(&vh[i], sizeof(NexVarHeader), 1, fp);
	}

    // 3. read timestamp data for the first variable
	// assuming the first variable is neuron or event
	// and it contains less than 10000 timestamps
	int timestamps[10000];
	if(vh[0].Type == 0 || vh[0].Type == 1){
		// seek to the start of data
		fseek(fp,  vh[0].DataOffset, SEEK_SET);
		// read the timestamps, 4 bytes per timestamp
		fread(timestamps, vh[0].Count*4, 1, fp);
	}
	fclose(fp);
  }
  -------------------------------------------------------


  2. NEX FILE HEADERS
  -------------------
  
  2.1. NexFileHeader contains the following information:
  (full definition of NexFileHeader is provided below)

	int  Version; // 100 to 104
	char Comment[256]; // file comment
	double Frequency;  // timestamp frequency - tics per second
	int  Beg;	// usually 0, minimum of all the timestamps in the file
	int  End;	// = maximum timestamp + 1
	int  NumVars;    // number of variables (and variable headers) in the file

    NexFileHeader is followed by variable headers.

  2.2. NexVarHeader contains the following information:
  (full definition of NexVarHeader is provided below)

  	int  Type; // 0-neuron, 1-event, 2-interval, 3-waveform, 
	           // 4-population vector, 5-continuous variable, 6 - marker
	char Name[64]; // variable name
	int  DataOffset; // where the data array for this variable 
					//is located in the file
	int  Count;  // number of events, intervals, waveforms, or fragments
	double XPos; // neuron only, electrode position in (0,100) range
	double YPos; //  neuron only, electrode position in (0,100) range
	double WFrequency; // waveforms and continuous vars only, 
					// sampling frequency
	double ADtoMV; // waveforms and continuous vars only, 
				//coeff. to convert from A/D values to milliVolts.
	int  NPointsWave; // waveforms and continuous vars. only, 
				// waveforms: number of points in each wave
				// continuous vars.: total number of a/d valus

	int	 NMarkers; // how many values are associated with each marker
	int	 MarkerLength; // how many characters are in each marker value


  3. HOW DATA ARE STORED
  ---------------------

	3.1. Neurons and events (VarHeader.Type = 0 or VarHeader.Type = 1).
	
	  The timestamps are stored as arrays of 4-byte integers.
		NexVarHeader vh;
		int timestamps[10000];
		// seek to the start of data
		fseek(fp,  vh.DataOffset, SEEK_SET);
		// read the timestamps, 4 bytes per timestamp
		fread(timestamps, vh.Count*4, 1, fp);

    3.2. Interval variables (VarHeader.Type = 2)

	  Interval beginnings and interval ends are stored as 2 
	  arrays of 4-byte integers.
		NexVarHeader vh;
		int beginnings[10000], ends[10000];
		// seek the start of data
		fseek(fp,  vh.DataOffset, SEEK_SET);
		// read interval beginnings and ends
		fread(starts, vh.Count*4, 1, fp);
		fread(ends, vh.Count*4, 1, fp);
		
     3.3  Waveform variables (VarHeader.Type = 3)

	  The data for a waveform variable are stored as:
	  - array of timestamps (4-byte integers)
	  - array of waveform values (2-byte integers, raw A/D values) 
	    for all the waveforms
		
		NexVarHeader vh;
		int timestamps[10000];
		short waveforms[32*10000];
		fseek(fp,  vh.DataOffset, SEEK_SET);
		fread(timestamps, vh.Count*4, 1, fp);
		fread(waveforms, vh.Count*vh.NPointsWave*2, 1, fp);

  		You also need to use the following fields in the variable header:
		vh.WFrequency = 25000; // this is a/d frequency of the waveform values, 25 kHz
		vh.ADtoMV = 1.; // this is a coefficient to convert a/d values to millivolts. 
		// 1. here means that the stored a/d values are in millivolts. 

      3.4  Continuously recorded variables (VarHeader.Type = 5)

      In general, a continuous variable may contain several fragments of data.
      Each fragment may be of different length. We don't store the timestamps
	  for all the a/d values since they would use too much space. Instead,
	  for each fragment we store the timestamp of the first a/d value in
	  the fragment and the index of the first data point in the fragment.
	  
	  Therefore, a continous variable contains the following 3 arrays:

      - array of all a/d values (stored as 2-byte integers) 
		NexVarHeader.NPointsWave field stores the number of a/d values

      - array of timestamps (each timestamp is for the beginning of the fragment;
        timestamps are in the same units as any other timestamps, i.e. usually
        in 25 usec units for Plexon data)
		NexVarHeader.Count field stores the number of fragments

      - array of indexes (each index is the position of the first data point
        of the fragment in the a/d array; index[0] is always 0, if index[1] = 200,
        it means that the second fragment is advalues[200], advalues[201], etc.)

		NexVarHeader vh;
		// assume that we have less than 100 fragments and less than 10000
		// a/d values
		int fragment_timestamps[100];
		int fragment_indexes[100];
		short advalues[10000];
		fseek(fp,  vh.DataOffset, SEEK_SET);
		fread(fragment_timestamps, vh.Count*4, 1, fp);
		fread(fragment_indexes, vh.Count*4, 1, fp);
		fread(advalues, vh.NPointsWave*2, 1, fp);

		You also need to use the following fields in the variable header:
		vh.WFrequency = 10; // this is a/d frequency, 10 data points per second
		vh.ADtoMV = 1; // this is a coefficient to convert a/d values to millivolts. 
		// 1. here means that the stored a/d values are in millivolts. 



      3.5 Markers (VarHeader.Type = 6)

       In general, a marker variable may contain several data fields, i.e.
       for each timestamp there are several values associated with this
       timestamp. The values of these data fields are stored as strings.

       If you are using Plexon system, there is only one data field for
       each timestamp (the strobed value) and the numerical strobed values
       are stored as strings. If you need the numerical strobed values,
       you'll have to convert strings to numbers.

		The data are stored in the following way:

		- array of timestamps
		for each field:
			- field names (each field name uses 64 bytes in the file)
			- array of field values (each value uses MarkerLength bytes in the file)
  
		When you read variable headers, check for Type = 6

	    Use then the following NexVarHeader fields:

			Name - variable name
			Count - number of timestamps
			NMarkers - number of data fields for each timestamp
			MarkerLength - length of each data field string 
			DataOffset - where the data starts in the file

		  Here is the code to read the marker variable:

  			// read the variable header 
			NexVarHeader sh;
			int* timestamps = 0;
			// assume that we have less than 16 fields
			char fieldnames[16][64];
			// assume that marker strings are shorter than 1024 characters
			char buf[1024];
			fread(&sh, sizeof(NexVarHeader), 1, fp);
			if(sh.Type == 6){
				timestamps = new int[sh.Count];
			}

			// seek to the variable data
			fseek(fp,  sh.DataOffset, SEEK_SET);

			// read the timestamps
			fread(timestamps, sh.Count*4, 1, fp);

			// read the field names and values
			for(i=0; i<sh.NMarkers; i++){
				// read the name of the data fields
				fread(fieldnames[i], 64, 1, fp);
				for(j=0; j<sh.Count; j++){
					// this is the i-th field for j-th timestamp
					fread(buf, sh.MarkerLength, 1, fp);
				}
			}

  
	4. Header Sizes
    ---------------

	Header sizes (to check integer sizes and structure alignment):

	sizeof(NexFileHeader) = 544
	sizeof(NexVarHeader) = 208

*/

#include "stdio.h"
#include "string.h"

// Nex file header structure
struct NexFileHeader {
	int  MagicNumber; // string NEX1
	int  Version; // 100
	char Comment[256];
	double Frequency;  // timestamped freq. - tics per second
	int  Beg;	// usually 0
	int  End;	// = maximum timestamp + 1
	int  NumVars;    // number of variables in the first batch
	int  NextFileHeader; // position of the next file header in the file
	//, not implemented yet
	char Padding[256]; // future expansion
};

// Nex variable header structure
struct NexVarHeader {
	int  Type; // 0 - neuron, 1 event, 2- interval, 3 - waveform, 4 - pop. vector, 5 - continuously recorded
	int  Version; // 100
	char Name[64]; // variable name
	int  DataOffset; // where the data array for this variable is located in the file
	int  Count;  // number of events, intervals, waveforms or weights
	int  WireNumber; // neuron only, not used now
	int  UnitNumber; // neuron only, not used now
	int  Gain;        // neuron only, not used now
	int  Filter;        // neuron only, not used now
	double XPos; // neuron only, electrode position in (0,100) range, used in 3D
	double YPos; //  neuron only, electrode position in (0,100) range, used in 3D
	double WFrequency; // waveform and continuous vars only, w/f sampling frequency
	double ADtoMV; // waveform continuous vars only, coeff. to convert from A/D values to Millivolts.
	int  NPointsWave; // waveform only, number of points in each wave
	int		NMarkers; // how many values are associated with each marker
	int		MarkerLength; // how many characters are in each marker value
	char	Padding[68];
};

int WriteNexFile(const char* fname);
int ReadNexFile(const char* fname);


int main()
{
	// write a file
	WriteNexFile("test.nex");
	// try to read a nex file
	ReadNexFile("test.nex");

	return 0;
}


int WriteNexFile(const char* fname)
{
	
	FILE* fp;

	fp = fopen("test.nex", "wb");
	if(fp == 0) return -1;

	// we'll make a nex file with 6 variables:
	// - 2 spike trains 
	// - 1 interval variable
	// - 1 continuous variable
	// - 1 waveform variable
	// - 1 marker variable

	// first, let's fill the spike trains
	// we will work with time resolution of 100 us or frequency 10000 Hz
	// neuron 1: spikes at 1, 2, and 3 sec
	int train1[3] = { 10000, 20000, 30000} ; // in time ticks: 1.*10000, etc.
	// neuron 2: spikes at 5, and 6 sec
	int train2[2] = { 50000, 60000} ; // in time ticks
	// interval variable with 2 intervals: (1 sec, 2 sec) and (4 sec, 6 sec)
	int intstart[2] = { 10000, 40000 };
	int intend[2]   = { 20000, 60000 };

	// continuous variable: 2 fragments, one starts at 1 sec, another at 7 sec
	int num_fragments = 2;
	int cont_ts[2] = { 10000, 70000 };
	// continuous variable: first fragment has 5 data points, second has 3 data points
	// what we need to store is the index of the first data point in each fragment
	// for the first fragment, the first data point is point number 0, 
	// for the second fragment, the first data point is point number 5 
	int cont_num_points[2] = { 0, 5 };
	// continuous variable: a/d values
	int num_data_points = 8;
	short cont_data_points[8] = { 0, 1, 2, 3, 4, 0, 4, 6 };

	// waveform variable: 3 waveforms, start at 0.5 sec, 1.5 sec and 2 sec
	int num_wf = 3;
	int wf_ts[3] = { 5000, 15000, 20000 };
	// waveform variable: a/d values
	int num_data_points_in_wave = 8;
	// waveform values
	short waveform_values1[8] = { 0,1,2,3,4,5,6,7};
	short waveform_values2[8] = { 7,6,5,4,3,2,1,0};
	short waveform_values3[8] = { 0,1,2,3,4,5,6,7};

	// marker variable: 3 markers, each has 2 fields
	// markers are at 1, 2, and 3 sec
	int marker_ts[3] = { 10000, 20000, 30000} ; // in time ticks: 1.*10000, etc.
	// names of marker fields;
	char marker_field_names[2][64];
	memset(marker_field_names[0], 0, 64);
	strcpy(marker_field_names[0], "field1");
	memset(marker_field_names[1], 0, 64);
	strcpy(marker_field_names[1], "field2");

	// values of the first marker field:
	char* marker_field_values1[] = {"field1value1", "field1value2", "field1value3"};
	// values of the second marker field:
	char* marker_field_values2[] = {"field2value1", "field2value2", "field2value3"};


	NexFileHeader fh;
	memset(&fh, 0, sizeof(NexFileHeader));

	char magic[] = "NEX1";
	fh.MagicNumber = *(int*)magic;
	fh.Version = 100;
	strcpy(fh.Comment, "test nex file");
	fh.Frequency = 10000.; // frequency is 10000 hz, or 100 us time resolution
	fh.Beg = 0;
	fh.End = 10000*10; // 10 seconds of data in the file
	fh.NumVars = 6;  // we have 6 variables in the file  
	fh.NextFileHeader = 0; 

	fwrite(&fh, sizeof(NexFileHeader), 1, fp);

	// calculate the data offset for the first variable
	// right after all the headers
	int offset = sizeof(NexFileHeader) + fh.NumVars*sizeof(NexVarHeader);

	NexVarHeader varheader;
	memset(&varheader, 0, sizeof(NexVarHeader));

	// neuron 1 header
	varheader.Type = 0; // neuron type
	varheader.Version = 100; 
	strcpy(varheader.Name, "neuron1"); // variable name
	varheader.DataOffset = offset;
	varheader.Count = 3;  // number of events
	fwrite(&varheader, sizeof(NexVarHeader), 1, fp);

	// adjust offset: add the bytes of the neuron1 timestamp array
	offset += varheader.Count*sizeof(int);
	
	// neuron 2 header
	memset(&varheader, 0, sizeof(NexVarHeader));
	varheader.Type = 0; // neuron type
	varheader.Version = 100; 
	strcpy(varheader.Name, "neuron2"); // variable name
	varheader.DataOffset = offset;
	varheader.Count = 2;  // number of fragments
	fwrite(&varheader, sizeof(NexVarHeader), 1, fp);

	// adjust offset: add the bytes of the neuron2 timestamp array
	offset += varheader.Count*sizeof(int);

	// interval header
	memset(&varheader, 0, sizeof(NexVarHeader));
	varheader.Type = 2; // interval type
	varheader.Version = 100; 
	strcpy(varheader.Name, "interval1"); // variable name
	varheader.DataOffset = offset;
	varheader.Count = 2;  // number of intervals
	fwrite(&varheader, sizeof(NexVarHeader), 1, fp);

	// adjust offset: add the bytes of the interval starts and ends
	offset += varheader.Count*sizeof(int)*2;

	// continuous variable header
	memset(&varheader, 0, sizeof(NexVarHeader));
	varheader.Type = 5; // continuous type
	varheader.Version = 100; 
	strcpy(varheader.Name, "continuous1"); // variable name
	varheader.DataOffset = offset;
	varheader.Count = 2;  // number of fragments
	varheader.NPointsWave = 8; // we have 8 a/d values in 2 fragments
	varheader.WFrequency = 10; // this is a/d frequency, 10 data points per second
	varheader.ADtoMV = 0.001; // this is a coefficient to convert a/d values ot millivolts. 
	// 0.001 here means that the stored a/d values are in volts. 

	fwrite(&varheader, sizeof(NexVarHeader), 1, fp);

	// adjust offset: add the bytes of fragment timestamps
	offset += varheader.Count*sizeof(int);
	// adjust offset: add the bytes of fragment number of values
	offset += varheader.Count*sizeof(int);
	// adjust offset: add the bytes of a/d values
	offset += varheader.NPointsWave*sizeof(short);

	// waveform variable header
	memset(&varheader, 0, sizeof(NexVarHeader));
	varheader.Type = 3; // waveform type
	varheader.Version = 100; 
	strcpy(varheader.Name, "wave1"); // variable name
	varheader.DataOffset = offset;
	varheader.Count = 3;  // number of waveforms
	varheader.NPointsWave = 8; // we have 8 a/d values in each waveform
	varheader.WFrequency = 1000; // this is a/d frequency, 1000 data points per second
	varheader.ADtoMV = 1.; // this is a coefficient to convert a/d values ot millivolts. 
	// 1 here means that the stored a/d values are in millivolts. 

	fwrite(&varheader, sizeof(NexVarHeader), 1, fp);

	// adjust offset: add the bytes of waveform timestamps
	offset += varheader.Count*sizeof(int);
	// adjust offset: add the bytes of a/d values
	offset += varheader.Count*varheader.NPointsWave*sizeof(short);

	// marker header
	memset(&varheader, 0, sizeof(NexVarHeader));
	varheader.Type = 6; // marker type
	varheader.Version = 100; 
	strcpy(varheader.Name, "marker1"); // variable name
	varheader.DataOffset = offset;
	varheader.Count = 3;  // number of markers (timestamps)
	varheader.NMarkers = 2; // number of marker fields for each timestamp
	varheader.MarkerLength = 16; // number of bytes ofr each field string
	fwrite(&varheader, sizeof(NexVarHeader), 1, fp);

	// adjust offset: add the bytes of timestamps
	offset += varheader.Count*sizeof(int);
	// adjust offset: add the bytes of field names (64 bytes per field) and field values
	offset += varheader.NMarkers*(64 + varheader.Count*varheader.MarkerLength);
	
	// now, write the data in the same order as we wrote the headers
	// first spike train
	fwrite(train1, 3*sizeof(int), 1, fp);
	// second spike train
	fwrite(train2, 2*sizeof(int), 1, fp);
	// interval variable -- for the intervals, first we write starts, then ends
	fwrite(intstart, 2*sizeof(int), 1, fp);
	fwrite(intend, 2*sizeof(int), 1, fp);

	// for the continuous variable, first we write fragment timestamps, 
	// then fragment counts, and finally, a/d values
	fwrite(cont_ts, num_fragments*sizeof(int), 1, fp);
	fwrite(cont_num_points, num_fragments*sizeof(int), 1, fp);
	fwrite(cont_data_points, num_data_points*sizeof(short), 1, fp);

	// for waveform variable, first we write w/f timestamps, 
	// then, a/d values
	fwrite(wf_ts, num_wf*sizeof(int), 1, fp);
	fwrite(waveform_values1, 8*sizeof(short), 1, fp);
	fwrite(waveform_values2, 8*sizeof(short), 1, fp);
	fwrite(waveform_values3, 8*sizeof(short), 1, fp);

	
	char buf[16];
	int i;
	// for the marker variable, first we write the timestamps, 
	fwrite(marker_ts, 3*sizeof(int), 1, fp);
	// then, for each field -- 
	// field name (we have to use exactly 64 bytes for each name)
	fwrite(marker_field_names[0], 64, 1, fp);
	// then, marker field values for the first field
	// we have to use exactly varheader.MarkerLength bytes for each value
	for(i=0; i<3; i++){
		strcpy(buf, marker_field_values1[i]);
		fwrite(buf, 16, 1, fp);
	}
	// field name for the second filed
	fwrite(marker_field_names[1], 64, 1, fp);
	// finally, marker field values for the second field
	// we have to use exactly varheader.MarkerLength bytes for each value
	for(i=0; i<3; i++){
		strcpy(buf, marker_field_values2[i]);
		fwrite(buf, 16, 1, fp);
	}

	printf("test.nex created!\n");
	fclose(fp);
	return 0;
}

int ReadNexFile(const char* fname)
{

	FILE* fp;

	fp = fopen(fname, "rb");
	if(fp == 0){
		printf("cannot open file %s\n", fname);
		return -1;
	}

	// 1. read the file header
	NexFileHeader fh;
	fread(&fh, sizeof(NexFileHeader), 1, fp);

	printf("nex file version: %d\n", fh.Version);
	printf("frequency: %f hz\n", fh.Frequency);
	printf("duration: %f sec\n", (double)fh.End/fh.Frequency);
	printf("number of variables: %d\n\n", fh.NumVars);
	printf("variables:\n\n", fh.NumVars);
	
	
	// 2. read the variable headers
	NexVarHeader varheaders[100];
	fread(varheaders, sizeof(NexVarHeader), fh.NumVars, fp);
	int i, j, k;
	for(i=0; i<fh.NumVars; i++){
		printf("%s type %d count %d\n", varheaders[i].Name, varheaders[i].Type, 
			varheaders[i].Count);
	}

	// 3. read the variable data
	// we'll read first two neurons or events only
	// assume that each neuron or event has less than 100000 timestamps
	int timestamps[100000];
	int buffer[100000];
	short sbuffer[100000];
	for(i=0; i<fh.NumVars; i++){
		if(varheaders[i].Type == 0 || varheaders[i].Type == 1){ // neuron or event
			// seek to the start of data
			fseek(fp,  varheaders[i].DataOffset, SEEK_SET);
			// read the timestamps
			fread(timestamps, varheaders[i].Count*4, 1, fp);
			// print first 5 timestamps as values in seconds
			printf("%s timestamps: ", varheaders[i].Name);
			for(j=0; j<varheaders[i].Count; j++){
				if(j >= 5) break;
				printf("%.6f ", (double)timestamps[j]/fh.Frequency);
			}
			printf("\n");
		}
		else if(varheaders[i].Type == 2){ // interval
			// seek to the start of data
			fseek(fp,  varheaders[i].DataOffset, SEEK_SET);
			// read the timestamps for interval starts
			fread(timestamps, varheaders[i].Count*4, 1, fp);
			// read the timestamps for interval ends
			fread(buffer, varheaders[i].Count*4, 1, fp);
			// print first 5 intervals as values in seconds
			printf("%s intervals: ", varheaders[i].Name);
			for(j=0; j<varheaders[i].Count; j++){
				if(j >= 5) break;
				printf("(%.6f, %.6f) ", (double)timestamps[j]/fh.Frequency,
					(double)buffer[j]/fh.Frequency);
			}
			printf("\n");
		}
		else if(varheaders[i].Type == 3){ // waveform variable
			// seek to the start of data
			fseek(fp,  varheaders[i].DataOffset, SEEK_SET);
			// read the timestamps for waveforms starts
			fread(timestamps, varheaders[i].Count*4, 1, fp);
			// read the a/d values of waveforms
			fread(sbuffer, varheaders[i].Count*varheaders[i].NPointsWave*2, 1, fp);
			printf("%s, %d waveforms\n", varheaders[i].Name, varheaders[i].Count);
			for(j=0; j<varheaders[i].Count; j++){
				if(j >= 5) break;
				printf("waveform %d starts at %.6f, values: ", j+1, 
					(double)timestamps[j]/fh.Frequency);
				int kstart = j*varheaders[i].NPointsWave;
				int kend = (j+1)*varheaders[i].NPointsWave;
				for(k = kstart; k< kend; k++){
					if(k - kstart > 3) break;
					printf("%.4f,", sbuffer[k]*varheaders[i].ADtoMV);
				}
				printf("\n");
			}
			printf("\n");

		}

		else if(varheaders[i].Type == 5){ // continuous variable
			// seek to the start of data
			fseek(fp,  varheaders[i].DataOffset, SEEK_SET);
			// read the timestamps for fragment starts
			fread(timestamps, varheaders[i].Count*4, 1, fp);
			// read the indexes for intervals
			fread(buffer, varheaders[i].Count*4, 1, fp);
			// read the a/d values
			fread(sbuffer, varheaders[i].NPointsWave*2, 1, fp);
			// print first 5 fragments
			printf("%s fragments:\n", varheaders[i].Name);
			for(j=0; j<varheaders[i].Count; j++){
				if(j >= 5) break;
				printf("fragment %d starts at %.6f, index %d; values: ", j+1, (double)timestamps[j]/fh.Frequency,
					buffer[j]);
				// print a few values in this fragment
				int kstart = buffer[j];
				int kend;
				if(j < varheaders[i].Count-1) kend = buffer[j+1];
				else kend = varheaders[i].NPointsWave;
				// print first 3 values of the fragment
				for(k = kstart; k< kend; k++){
					if(k - kstart > 3) break;
					printf("%.4f,", sbuffer[k]*varheaders[i].ADtoMV);
				}
				printf("\n");
			}
			printf("\n");
		}
		else if(varheaders[i].Type == 6){ // marker variable
			// seek to the start of data
			fseek(fp,  varheaders[i].DataOffset, SEEK_SET);
			// read the timestamps 
			fread(timestamps, varheaders[i].Count*4, 1, fp);
			// print first 5 timestamps as values in seconds
			printf("%s timestamps: ", varheaders[i].Name);
			for(j=0; j<varheaders[i].Count; j++){
				if(j >= 5) break;
				printf("%.6f ", (double)timestamps[j]/fh.Frequency);
			}
			printf("\n");
			printf("%s fields:\n", varheaders[i].Name);
			// read each field information
			for(j=0; j<varheaders[i].NMarkers; j++){
				char buf[1024];
				// read the name of the field
				fread(buf, 64, 1, fp);
				printf("field %d name=%s; field values: ", j+1, buf);
				// read the field values
				for(k=0; k<varheaders[i].Count; k++){
					if(k > 5) break;
					fread(buf, varheaders[i].MarkerLength, 1, fp);
					printf("%s, ", buf);
				}

				printf("\n");
			}
			printf("\n");
		}
	}
	printf("\n");
	fclose(fp);
	return 0;
}

