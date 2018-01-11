# Activity-Recognition
SVM classification for the MSR daily action dataset.

1) Code for RAD representation - RAD.cpp
2) Code for HJPD representation - HJPD.cpp
3) Code for HOD representation - HOD.cpp
4) Output for RAD - rad_d1 and rad_d1.t
5) Output for HJPD - hjpd_d1 and hjpd_d1.t
5) Output for HOD - hod_d1 and hod_d1.t

**Copy the 3 cpp files into a new folder in 'Home' called 'codes'
**We read the full dataset for all codes.
The dataset can be found from the course website, make sure to use 'dataset_full.tar.gz'
This is to be extracted into a new folder also in 'Home', called 'dataset_full'  
All codes create a file from where the file name is taken one at at time, called 'filenames.txt'.

RAD:
Complie code by following these steps-
Open a new terminal.
cd codes
g++ RAD.cpp
./a.out
The number of bins selected is 5 for both distance and angle.
The bin size is calculated by dividing the range of values by the number of bins.
Histograms X axis values start from 0, as we subtracted the min value from each of the value. This makes it easier to view. 
To run the Train and Test dataset, lines/parameters to be changed are shown clearly in the code.


HJPD:
Complie code by following these steps-
Open a new terminal.
cd codes
g++ HJPD.cpp
./a.out
The number of bins selected is 10.
The bin size is calculated by dividing the range of values by the number of bins.
Histograms X axis values start from 0, as we subtracted the min value from each of the value. This makes it easier to view. 
To run the Train and Test dataset, lines/parameters to be changed are shown clearly in the code.


HOD:
Complie code by following these steps-
Open a new terminal.
cd codes
g++ HJPD.cpp
./a.out
The 360 angle is divided into 8 parts of 45 degrees each. Further calculations are done accordingly.
Temporal pyramids are built as shown in class to the third level, so 7 pyramids each.
To run the Train and Test dataset, lines/parameters to be changed are shown clearly in the code.
