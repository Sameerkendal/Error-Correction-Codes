README:


The .zip file contains the following files :

1. RS_decoder.m : It's the main file to load corrupted mp3 file blockwise and decode it using berlekamp_welch_decoder function blockwise. It also rearrange the decoded blocks into mp3 file.

2.berlekamp_welch_decoder : It's the main function which takes corrupted signal as input and decode it back considering the RS coding. It implements berlekamp welch decoding algorithm to recover the actual message. It uses gaussian elimination method to solve the equation and has some helper functions to compute the polynomial division.
All edge cases like 'A' being non-full rank matrix and padding zeros has been taken into  consideration.

3. Example.m : This is the example file which shows our implementation on lower scale and can be tuned accordingly. The messages and errors are hard-coded to check each step and satisfy maths carefully.

4.music8_corrupted.mp3 : This is corrupted mp3 file provided to us as an input.

5. Group8_decoded.mp3 : This is corrected mp3 file generated ater running RS_decoder.m file to decode the input file.

