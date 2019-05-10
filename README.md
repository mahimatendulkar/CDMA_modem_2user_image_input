# CDMA_modem_2user_image_input
Code for CDMA transmission and reception for 2 users (image transmission)



ALGORITHM USED
1.	Took different images as input 1 and input 2.
2.	These images where then resized and reshaped suitably.
3.	Since both of them produce different length arrays, Zero padding operation is performed.
4.	Now the user1’s and user2’s  input is converted to Bipolar NRZ
5.	Assumed the carrier frequency as 1, Energy per bit as 2, and time per bit of message sequence as 1.
6.	Designed the CDMA transmitter for user1 and user2.
7.	BPSK modulation was implemented for both.
8.	Generated PN sequences and multiplied with the modulated output.
9.	Composite signal was obtained by adding outputs of both users from above step.
10.	This signal was then passed through AWGN channel, with SNR(dB)=50.
11.	Signal received was then recovered and correlated with the code.
12.	BPSK demodulation was implemented and image was reconstructed.
