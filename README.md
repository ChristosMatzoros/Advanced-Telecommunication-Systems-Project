# ECE511-Advanced-Telecommunication-Systems
# Telecommunication system simulator in Matlab

The project is composed of four sub-assignments that use the Monte-Carlo method to simulate a telecommunication system (in Matlab).

In the first assignment, I implemented from scratch a discrete digital transmission AWGN channel with m-PSK modulation(BPSK/QPSK/8-PSK/16-PSK), with a changeable SNR value, size of oversampling and other parameters. It includes a BER evaluation for different SNR values.

In the second assignment, I appended the capability of simulating a fading channel(y=xh+w) with a varying channel characteristics (h). To obtain the correct values, the receiver implements equalization/channel inversion. The receiver can either know or not the value of h. In addition to that, I implemented the Maximal Ratio Combining(MRC) diversity technique for 2 diversity branches. A BER evaluation with different values of SNR is also included for every technique.

I continued by simulating a 2x2 MIMO fading channel for a varying channel. I implemented the spatial multiplexing technique, with a receiver that uses Least Squares channel inversion (2x2 MIMO). As a next step, I implemented a 16-QAM Modulation for the simulation. An evaluation was made to specify the best combination that maximizes the goodput.

Finally, on the last assignment, I simulated a simple version of an OFDM system (ΙΕΕΕ 802.11a) with a changeable cyclic prefix for a fading channel as well as for a steady LTI channel.
