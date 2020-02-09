%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Advanced  Topics        %  
%  In Telecommunication Systems  %
%           Homework 4           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matzoros Christos Konstantinos % 
%            AEM:2169            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%packages
pkg load signal;
pkg load communications;
close all;

%%%%%%%%%%%%%%%%%%%%%%%
% FOR HOMEWORKS 1 - 3 %
%%%%%%%%%%%%%%%%%%%%%%%
%Instructions in order to run the program correctly for every problem the homework 2:
spect_flag = 0;
%Uncomment the following line in order to compute and present the spectrum for the x,y of the y=x*h+w
%spect_flag = 1;

c_flag=0;
c_solution=0;
#Uncomment the following 3 lines in order to run the question c: 
%c_flag = 1;                            %flag to run the code for the 3rd question
%C = 10;                                 %C(constant): for the h'=h+C in the receiver 
%c_solution = 1;                        %uncomment this line in order to obtain the solution of problem c

d_flag=0;
%Uncomment the following 4 lines in order to run the question d:
%L=2;                                    %L = number of antennas on the receiver's side
%d_flag=1;
%spect_flag=0;
%c_flag = 0;

%Instructions in order to run the program correctly for every problem the homework 3:
MIMO = 0;
#Uncomment the following 5 lines in order to run the question a(Use of MIMO):
%L=2;                                    %L = 2 for MIMO(2x2)
%d_flag=1;
%spect_flag=0;
%c_flag = 0;
%MIMO = 1;

QAM_16 = 0;
%Uncomment the following 2 lines in order to run the question b(16-QAM):
#QAM_16 = 1;
#M_QAM = 16;                     % 16 QAM

GOODPUT_CALC = 0;
%Uncomment the following line to calculate goodput
%GOODPUT_CALC = 1;

%%%%%%%%%%%%%%%%%%
% FOR HOMEWORK 4 %
%%%%%%%%%%%%%%%%%%
%Instructions in order to run the program correctly for every problem the homework 4:
OFDM = 0;
#Uncomment the following lines in order to run the question a and b:
OFDM = 1;
NUM_OF_SUBCARRIERS = 64;
CP_SIZE = 16;

spect_flag = 0;
%Uncomment the following line in order to compute and present the spectrum
%spect_flag = 1;

LTI_OFDM = 0;
%Uncomment the following line in order to run question d (LTI system)
%LTI_OFDM = 1;

QAM_16 = 0;
%Uncomment the following 2 lines in order to run the questions with (16-QAM):
#QAM_16 = 1;
#M_QAM = 16; 


%Constants 
Rb = 10000;
m = 4;                              %m-PSK m = 2,4,8,16
N= 13312;                       
os = 1;                             %oversampling    
Eb = 1;
SNRdB = [0 2 4 6 8 10 12 14];
#SNRdB = [4];
Tchannel = 1;                     %Channel time period in seconds
Tpacket = 1; 

%other variables
Tb = 1/Rb;
num_of_symbols = floor(N/log2(m));
im_j = sqrt(-1);
errors_array = [];
BER = [];
BERavg = [];
SNRSlin = [];
mean=0;                                             %mean of the noise distribution
Ts = log2(m) * Tb;
Fs = 1/Ts;
mn = [];


function y = myqammodulator (x, m, Eb) 
  b = -2 .* mod (x, (sqrt(m))) + sqrt(m) - 1;
  a = 2 .* floor (x ./ (sqrt(m))) - sqrt(m) + 1;
  y = (sqrt(Eb))*(a + i.*b);

endfunction

function z = myqamdemodulator(y,m,Eb)
    c = sqrt(m);
    x = myqammodulator(0:(m-1),m,Eb);
    x = reshape(x,1,m);
    for k = 1:length(y)
        [n z(k)] = min(abs(y(k) - x));
        z(k) = z(k) - 1;
    end
endfunction

%this function is used in order to create and present the PSD/spectrum of the OFDM signal
function k = OFDM_PSD(nBitPerSymbol,nSymbol,nBit,input_bits,m,Fs,snrdb) 
    nFFTSize = 64;
    cp = nFFTSize/4;       %(25%)
    nBitPerSymbol = 52;
    
    sFilt = pskmod(input_bits,m);
    rest_bits = mod(size(sFilt,2),nBitPerSymbol);
    i=size(sFilt,2);

    s = i;
    while(i>(s-rest_bits))
        sFilt(i) = [];
        i=i-1;
    end        
    num_of_symbols =  size(sFilt,2);
            
    nBit  = num_of_symbols;
    nSymbol = ceil(nBit/nBitPerSymbol);
    subcarrierIndex = [-ceil(nBitPerSymbol/2):-1 1:ceil(nBitPerSymbol/2)];
    nSymbol = ceil(nBit/nBitPerSymbol); #number of OFDM symbols
    ip = sFilt;    
    ipMod = [ip zeros(1,nBitPerSymbol*nSymbol-nBit);];
    ipMod = reshape(ipMod,nSymbol,nBitPerSymbol);

    st = []; % empty vector
    for i = 1:nSymbol
        inputiFFT = zeros(1,nFFTSize);
        % assigning bits a1 to a52 to subcarriers [-26 to -1, 1 to 26]
        inputiFFT(subcarrierIndex+nFFTSize/2+1) = ipMod(i,:);
        
        %  shift subcarriers at indices [-26 to -1] to fft input indices [38 to 63]
        inputiFFT = fftshift(inputiFFT);
        outputiFFT = ifft(inputiFFT,nFFTSize);
        % adding cyclic prefix of 16 samples 
        
        cycle_prefix = outputiFFT(nFFTSize-cp+1 : nFFTSize);
        outputiFFT_with_CP = [cycle_prefix outputiFFT];   
        st = [st outputiFFT_with_CP];
    end
    
    st = awgn(st,snrdb,'measured');
    
    figure(4)
    temp_st = st;
    fsMHz = Fs
    DFT_points = ceil(2*Fs);
    [Pxx,W] = pwelch(temp_st,[],[],DFT_points,fsMHz);    
    plot([-(DFT_points/2):((DFT_points/2 -1))]*fsMHz/DFT_points,10*log10(fftshift(Pxx)),"linewidth", 2);
    xlabel('frequency, MHz')
    ylabel('power spectral density')
    title('Transmit spectrum OFDM');
    xlim([-fsMHz/2 fsMHz/2])
    
    figure(1)
    
endfunction


screen_size = get(0, "screensize");
s_height = screen_size(4);
s_width = screen_size(3);

if(QAM_16 ==1)
    m = M_QAM;
end
if(m>4)
    GOODPUT_CALC = 0;
end

goodput = [];
theoretical_BER = [];
#this loop is used in order to make a monte carlo simulation with different values in SNR
snr_counter=1;
disp ("Press the space button or click on figure 2 in order to continue to the next value of SNR");
while(snr_counter<=size(SNRdB,2))
    num_of_symbols = floor(N/log2(m));
    SNRlin = 10^(SNRdB(snr_counter)/10);            %find SNRlin in order from the current SNR in dB
    No = (2*Eb)/SNRlin;                             %No = Noise
    variance = No/2;                                %variance(power) of noise
    sigma = sqrt(variance);                         %standard deviation
    th_ber =  (1/2) *(erfc(sqrt(SNRlin/2)));
    theoretical_BER = [theoretical_BER th_ber]; 
    
    %%%%%%%%%%%%%%%%%%%%
    %  Input Creation  %
    %%%%%%%%%%%%%%%%%%%%
    input_bits = randi([0 1],1,N);                  %create random input bits of ones and zeros
    bit_time = linspace(Tb,N*Tb,N);
    
    figure(1)
    set(figure(1), 'Position', [0 (s_height/2) (1.2*s_width/2) (2*s_height/5)])
    
    subplot (2, 3, 1);
    stem(bit_time,input_bits,"filled");
    title ("Input bits");

    %%%%%%%%%%%%%%%%%%%%%%
    %%  Symbol Encoder  %%
    %%%%%%%%%%%%%%%%%%%%%%
    
    if(QAM_16 == 0)
        %discard the excessive bits in order to create the symbols for the current m-PSK
        rest_bits = mod(N,log2(m));                             
        k = N-mod(N,log2(m));
        i=N;
        while(i>(N-mod(N,log2(m))))
            input_bits(i) = [];
            i=i-1;
        end
        N = size(input_bits,2);
        
        symbols = [];
        i=1;
        while(i<=N)
            bit_string = num2str([]);
            for k=0:log2(m)-1
                bit_string = strcat(bit_string,num2str(input_bits(i+k)));
            end
            l = bin2dec(bit_string);
            symbols = [symbols (sqrt(log2(m)*Eb) * exp(im_j*l*pi/(m/2) + im_j*(pi/m)))];    #symbol creation
            i=i+log2(m);
        end
        
           
        if(MIMO == 1)                       %In order for the symbols to be multiple of L in order to implement LxL MIMO
            rest_bits = mod(num_of_symbols,L);
            i=num_of_symbols;
            while(i>(num_of_symbols-rest_bits))
                symbols(i) = [];
                i=i-1;
            end
            
            num_of_symbols =  size(symbols,2);    
        end
        
    elseif(QAM_16 == 1)
        k_QAM = log2(M_QAM);                % Number of bits per symbol
        N_QAM = N;                          % Number of bits to process
        
        data = reshape(input_bits,size(input_bits,2)/k_QAM,k_QAM);   % Reshape data into binary k-tuples, k = log2(M)   
        symbs = [];
    
        number_of_decs = size(data,1);
        i=1;
        while(i<=number_of_decs)
            temp_tuple = data(i,:);
            
            bit_string = num2str([]);
            for j=1:log2(M_QAM)
                bit_string = strcat(bit_string,num2str(temp_tuple(j)));
            end
            l = bin2dec(bit_string);

            symbs = [symbs l];
            i=i+1;
        end

        symbols = myqammodulator(symbs,M_QAM,Eb); % Gray coding, phase offset = 0
        num_of_symbols = size(symbols,2);
        
        if(MIMO == 1)                       %In order for the symbols to be multiple of L in order to implement LxL MIMO
            rest_bits = mod(num_of_symbols,L);
            i=num_of_symbols;
            while(i>(num_of_symbols-rest_bits))
                symbols(i) = [];
                i=i-1;
            end
            
            num_of_symbols =  size(symbols,2);    
        end
        
    end
    
    re_s = real(symbols);
    im_s = imag(symbols);
    time_start = linspace(Ts,num_of_symbols*Ts,num_of_symbols);

    s1 = size(time_start);
    s2 = size(re_s);
    subplot (2, 3, 2);
    stem(time_start,re_s,"filled");
    title ("Real part of symbol encoder output");

    %%%%%%%%%%%%%%%%%
    %%  TX Filter  %%
    %%%%%%%%%%%%%%%%%
    re_sFilt = [];
    im_sFilt = [];
    filt_arr = (sqrt(Eb/os))*ones(1,os);
    for i=1:num_of_symbols
        re_value = re_s(i);
        im_value = im_s(i);
        re_sFilt =[re_sFilt (conv(re_value,filt_arr))];
        im_sFilt =[im_sFilt (conv(im_value,filt_arr))];
    end

    time = linspace(Ts/os,num_of_symbols*Ts,num_of_symbols*os);
    
    subplot (2, 3, 3);
    stem(time,re_sFilt,"filled");
    title ("Real part of TX filter output");
    
    %%%%%%%%%%%%%%%%%% For Homework 4 %%%%%%%%%%%%%%%%%%
    if(OFDM == 1)
        nFFTSize = NUM_OF_SUBCARRIERS;
        cp = CP_SIZE;       # (25%)
        nBitPerSymbol = 52;
        
        sFilt = re_sFilt + im_sFilt * im_j;
        rest_bits = mod(size(sFilt,2),nBitPerSymbol);
        i=size(sFilt,2);
        
        s = i;
        while(i>(s-rest_bits))
            sFilt(i) = [];
            i=i-1;
        end        
        num_of_symbols =  size(sFilt,2);
                
        nBit  = num_of_symbols;
        nSymbol = ceil(nBit/nBitPerSymbol);
        subcarrierIndex = [-ceil(nBitPerSymbol/2):-1 1:ceil(nBitPerSymbol/2)];
        nSymbol = ceil(nBit/nBitPerSymbol); #number of OFDM symbols
        
        ip = sFilt;
        ipMod = [ip zeros(1,nBitPerSymbol*nSymbol-nBit);];
        ipMod = reshape(ipMod,nSymbol,nBitPerSymbol);
         
        %%%%% For LTI channel
        if(LTI_OFDM==1)
            h= [0.9+0.9*im_j, 0.6+0.6*im_j,0.3+0.3*im_j];
            delay_spread = 3;
            cp = delay_spread;
            h_coef = [ h(1)'/abs(h(1)*h(1)')  h(2)'/abs(h(2)*h(2)')  h(3)'/abs(h(3)*h(3)')];                    
        end

        res = [];       
        for i = 1:nSymbol
            inputiFFT = zeros(1,nFFTSize);
            % assigning bits to subcarriers 
            inputiFFT(subcarrierIndex+nFFTSize/2+1) = ipMod(i,:); 
            %  shift subcarriers
            inputiFFT = fftshift(inputiFFT);           
            outputiFFT = ifft(inputiFFT,nFFTSize);        
            % adding cyclic prefix of 16 samples 
            cycle_prefix = outputiFFT(nFFTSize-cp+1 : nFFTSize);         
            outputiFFT_with_CP = [cycle_prefix outputiFFT];
            
            if(LTI_OFDM == 1)
                outputiFFT_with_CP = conv(outputiFFT_with_CP,h);                
                outputiFFT_with_CP = outputiFFT_with_CP(1:nFFTSize+cp);
            end
            
            res = [res outputiFFT_with_CP];
        end
        
        
        re_sFilt = real(res);
        im_sFilt = imag(res);
        num_of_symbols = size(re_sFilt,2);
        time = linspace(Ts/os,num_of_symbols*Ts,num_of_symbols*os);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      
    %%%%%%%%%%%%%%%%%%%%
    %%  Insert Noise  %%
    %%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%% For Homework 4 %%%%%%%%%%%%%%%%%%
    if(OFDM == 1)
        Noise_re = [];
        Noise_re = [Noise_re (sigma * randn(size(re_sFilt,2), 1))]+mean;
        Noise_re = reshape(Noise_re,1,size(re_sFilt,2));

        Noise_im = [];
        Noise_im = [Noise_im (sigma * randn(size(re_sFilt,2), 1))]+mean;
        Noise_im = reshape(Noise_im,1,size(re_sFilt,2));
        
       
        re_y = awgn(re_sFilt,SNRdB(snr_counter),'measured');
        im_y = awgn(im_sFilt,SNRdB(snr_counter),'measured');
        
        inputbits = re_y + im_y*im_j;
        
        if(spect_flag == 1)
            OFDM_PSD(nBitPerSymbol,nSymbol,nBit,input_bits,m,Fs,SNRdB(snr_counter));
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    elseif(MIMO == 1)
        l_counter = 1;
        noise_array_re=[];
        noise_array_im=[];
        while(l_counter<=L)
            Noise_re = [];
            Noise_re = [Noise_re (sigma * randn(size(re_sFilt,2), 1))]+mean;
            Noise_re = reshape(Noise_re,1,size(re_sFilt,2));
            
            Noise_im = [];
            Noise_im = [Noise_im (sigma * randn(size(re_sFilt,2), 1))]+mean;
            Noise_im = reshape(Noise_im,1,size(re_sFilt,2));
            
            noise_array_re = [noise_array_re ; Noise_re];
            noise_array_im = [noise_array_im ; Noise_im];
            l_counter=l_counter+1;
        end
      
        h_ant = [];
        l_counter=1;
        y_array=[];
        
        while(l_counter<=L)
            internal_counter = 0;
            
            while(internal_counter < L)
                timer = 0;
                counter = 0;
                stop_time = time(size(time,2));
                num_of_period = -1;
                h = [];
                l = 0;
                while(timer<stop_time/2)
                    counter=counter+1;    
                    timer=time(counter);
                    cur_num_of_period = ceil(timer/Tchannel);
                    if(cur_num_of_period>num_of_period)
                        num_of_period = cur_num_of_period;
                        l = l+1;
                        X = sqrt(1/sqrt(2))*randn(1) + mean;
                        Y = sqrt(1/sqrt(2))*randn(1)+ mean;
                        h_new = X + Y * im_j;
                    end
                    h =[h h_new]; 
                end
                h_ant=[h_ant;h];
            
                
                internal_counter = internal_counter + 1;
            end
            
            count = 1;
            y_eq = [];
            no = norm(h);
            y= [];
            
            index = (l_counter*L)-L+1;
            h_counter = 1;
            while(count <= size(re_sFilt,2))
                x1 =  re_sFilt(count)+im_sFilt(count)*im_j;
                x2 =  re_sFilt(count+1)+im_sFilt(count+1)*im_j;
                noise =  noise_array_re(l_counter,count)+ noise_array_im(l_counter,count)*im_j;
                y = [y (h_ant(index,h_counter)*x1+h_ant(index+1,h_counter)*x2+noise)];      
                h_counter = h_counter+1;
                count = count +2;
            end
        
            y_array = [y_array;y];
            l_counter = l_counter + 1;     
        end
  
    elseif((MIMO == 0) && (d_flag==0))   %for questions from a to c of homework 2
        Noise_re = [];
        Noise_re = [Noise_re (sigma * randn(size(re_sFilt,2), 1))]+mean;
        Noise_re = reshape(Noise_re,1,size(re_sFilt,2));

        Noise_im = [];
        Noise_im = [Noise_im (sigma * randn(size(re_sFilt,2), 1))]+mean;
        Noise_im = reshape(Noise_im,1,size(re_sFilt,2));

        timer = 0;
        counter = 0;
        stop_time = time(size(time,2));
        num_of_period = -1;
        h = [];
        
        l = 0;
        while(timer<stop_time)
            counter=counter+1;    
            timer=time(counter);
            cur_num_of_period = ceil(timer/Tchannel);
            if(cur_num_of_period>num_of_period)
                num_of_period = cur_num_of_period;
                l = l+1;
                X = sqrt(1/sqrt(2))*randn(1) + mean;
                Y = sqrt(1/sqrt(2))*randn(1)+ mean;
                h_new = X + Y * im_j;
            end
            h =[h h_new]; 
        end
        
        count = 1;
        y_eq = [];
        no = norm(h);
        y= [];
        while(count <= size(re_sFilt,2))
            x =  re_sFilt(count)+im_sFilt(count)*im_j;
            noise = Noise_re(count)+Noise_im(count)*im_j;
            y = [y (h(count)*x+noise)];
            count = count +1;
        end
        re_y = real(y);
        im_y = imag(y);
              
        if(spect_flag==1)
            figure(4);
            Rxx = xcorr(re_sFilt+im_sFilt*im_j,'biased');
            Rxxdft = abs(fft(Rxx));
            Rxxdft = abs(fftshift(fft(Rxx)));
            freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));
            plot(freq,Rxxdft);
            title ("Spectrum plot of x");
            xlabel("Frequency");
            ylabel("Power");
            
            figure(5);
            Rxx = xcorr(re_y+im_y *im_j,'biased');
            Rxxdft = abs(fft(Rxx));
            Rxxdft = abs(fftshift(fft(Rxx)));
            freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));
            plot(freq,Rxxdft);
            title ("Spectrum plot of y");
            xlabel("Frequency");
            ylabel("Power");
   
            figure(1)
        end
        
    elseif((MIMO == 0) && d_flag==1) %for the question D of homework 3
        
        l_counter = 1;
        noise_array_re=[];
        noise_array_im=[];
        while(l_counter<=L)
            Noise_re = [];
            Noise_re = [Noise_re (sigma * randn(size(re_sFilt,2), 1))]+mean;
            Noise_re = reshape(Noise_re,1,size(re_sFilt,2));
            
            Noise_im = [];
            Noise_im = [Noise_im (sigma * randn(size(re_sFilt,2), 1))]+mean;
            Noise_im = reshape(Noise_im,1,size(re_sFilt,2));
            
            noise_array_re = [noise_array_re ; Noise_re];
            noise_array_im = [noise_array_im ; Noise_im];
            l_counter=l_counter+1;
        end
      
        h_ant = [];
        l_counter=1;
        y_array=[];
        while(l_counter<=L)        
            timer = 0;
            counter = 0;
            stop_time = time(size(time,2));
            num_of_period = -1;
            h = [];
            l = 0;
            while(timer<stop_time)
                counter=counter+1;    
                timer=time(counter);
                cur_num_of_period = ceil(timer/Tchannel);
                if(cur_num_of_period>num_of_period)
                    num_of_period = cur_num_of_period;
                    l = l+1;
                    X = sqrt(1/sqrt(2))*randn(1) + mean;
                    Y = sqrt(1/sqrt(2))*randn(1)+ mean;
                    h_new = X + Y * im_j;
                end
                h =[h h_new]; 
            end
            h_ant=[h_ant;h];
       
            count = 1;
            y_eq = [];
            no = norm(h);
            y= [];
            while(count <= size(re_sFilt,2))
                x =  re_sFilt(count)+im_sFilt(count)*im_j;
                noise =  noise_array_re(l_counter,count)+ noise_array_im(l_counter,count)*im_j;
                y = [y (h(count)*x+noise)];
                count = count +1;
            end
        
            y_array = [y_array;y];
            l_counter = l_counter + 1;     
        end    
    end 
    
    %%%%%%%%%%%%%%%%%% For Homework 4 %%%%%%%%%%%%%%%%%%
    if(OFDM == 1)
        res = re_y+im_y*im_j;
        without_cp = [];
        for j = 1:nSymbol
            start_pos = ((j-1)*(nFFTSize+cp))+cp+1;
                
            
            if(LTI_OFDM == 1)          %for LTI channel
                y = res(start_pos-cp : start_pos + nFFTSize-1);
                
                # Decision-Feedback Equalizer
                z(1)= y(1)*h_coef(1);
                z(2)= (y(2)-z(1)*h(2))*h_coef(1);
                for l=3:nFFTSize+cp
                    z(l)=(y(l)-z(l-1)*h(2)-z(l-2)*h(3))*h_coef(1);
                endfor
                
                temp = z(cp+1:nFFTSize+cp);
                tempFFT = fft(temp,nFFTSize);
                tempFFT = fftshift(tempFFT);          
                pos = (nFFTSize-nBitPerSymbol)/2;
                temp_array = tempFFT(pos+1:(nFFTSize-pos+1));
                temp_array((nBitPerSymbol/2)+1) = [];
                without_cp = [without_cp; temp_array];
      
            else  
                temp = res(start_pos : start_pos + nFFTSize-1);
                tempFFT = fft(temp,nFFTSize); 
                tempFFT = fftshift(tempFFT);
                temp = (nFFTSize-nBitPerSymbol)/2;
                temp_array = tempFFT(temp+1:(nFFTSize-temp+1));
                temp_array((nBitPerSymbol/2)+1) = [];
                without_cp = [without_cp; temp_array];    
            end       
            
        end
        
        without_cp = reshape(without_cp,1,nBit);
        re_y = real(without_cp); 
        im_y = imag(without_cp); 
        
        num_of_symbols = size(re_y,2);
        time = linspace(Ts/os,num_of_symbols*Ts,num_of_symbols*os);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(d_flag==0)
        figure(1);
        subplot (2, 3, 4);
        stem(time,re_y,"filled");
        title ("Real part of channel output after the insertion of Noise");
    end

    
    %%%%%%%%%%%%%%%%%
    %%  RX FILTER  %%
    %%%%%%%%%%%%%%%%%
    
    if((d_flag==0) && (MIMO == 0))
        rx_conv_filter = (sqrt(Eb/os))*ones(1,os);  
        re_RX_output = conv(re_y,rx_conv_filter);
        im_RX_output = conv(im_y,rx_conv_filter);
        
        time = [time (linspace(num_of_symbols*Ts+Ts/os,num_of_symbols*Ts+(os-1)*Ts/os, os-1))];
        subplot (2, 3, 5);
        stem(time,re_RX_output,"filled");
        title ("Real part of RX filter output");
    
    else    %for question D of homework 2 and for MIMO of homework 3
        re_RX_output = [];
        im_RX_output = [];
        l_counter=1;
        while(l_counter<=L)
            rx_conv_filter = (sqrt(Eb/os))*ones(1,os);
            re_RX_output =[re_RX_output ;conv(real(y_array(l_counter,:)) ,rx_conv_filter)];
            im_RX_output =[im_RX_output ;conv(imag(y_array(l_counter,:)) ,rx_conv_filter)];
            l_counter = l_counter+1;
        end
    end        
    
    %%%%%%%%%%%%%%%%
    %%  SAMPLING  %%
    %%%%%%%%%%%%%%%%
    
    if(OFDM == 1)
        sampling_counter = os;
        real_sampling_output = [];
        imaginary_sampling_output = [];
        time_sampling = [];

        while(sampling_counter <= size(re_RX_output,2)-(floor(os/2)))
            real_sampling_output = [real_sampling_output re_RX_output(sampling_counter)];
            imaginary_sampling_output = [imaginary_sampling_output im_RX_output(sampling_counter)];
            sampling_counter = sampling_counter + os;
        end
        
        sampling_output = real_sampling_output+imaginary_sampling_output*im_j;
        time_start = linspace(Ts,size(sampling_output,2)*Ts,size(sampling_output,2));
     
    elseif((d_flag == 0) && (MIMO == 0))
        sampling_counter = os;
        real_sampling_output = [];
        imaginary_sampling_output = [];
        time_sampling = [];

        while(sampling_counter <= size(re_RX_output,2)-(floor(os/2)))
            real_sampling_output = [real_sampling_output re_RX_output(sampling_counter)];
            imaginary_sampling_output = [imaginary_sampling_output im_RX_output(sampling_counter)];
            sampling_counter = sampling_counter + os;
        end
       
        h_array = [];
        count = 1;
        index = 1;
        
        while(count <= size(real_sampling_output,2))
            h_array = [h_array h(index)];             
            count = count + 1;
            index = index + os;
        end
        
        temp_h = h_array;
        if(c_flag==1)  
            h_array = h_array+C;                                 %add the constant C that simulate the noise of a given channel h on the receiver (for question c)            
            %We know that our channel h is gaussian with mean=0. A new h' = h+C is gine to the receiver. We need to estimate the value C 
            %in order for the receiver to obtain the real value of h. So we have:
            %E[h'] = E[h+C] . C is a constant so,  E[h'] = E[h]+C. We know for our channel h, that mean=0 so E[h] = 0.
            %Finaly we get  E[h'] = C = mean of values of h'. In order to obtain an estimation for h, we substract the value C for every value of h'. 
            
            if(c_solution==1)
                mean_of_y = sum(real(h_array))/length(h_array);
                guess_of_C = mean_of_y;
                h_array = h_array - guess_of_C;
            end 
        end
        
        n_h_squared = [];
        y_eq = [];
        count = 1;
        X = real_sampling_output + imaginary_sampling_output * im_j;
        while(count <= size(real_sampling_output,2))
            n_h_squared =[n_h_squared (norm(temp_h(count)))^2];
            y_eq = [y_eq (X(count)*conj(h_array(count))/(norm(h_array(count))^2))];  % Y = X+(h*)/|h|^2 * n 
            count = count + 1;
        end
        
        real_sampling_output = real(y_eq);
        imaginary_sampling_output = imag(y_eq);
   
   
    elseif((d_flag == 1) && (MIMO == 0))           %for the question D of homework 2
        l_counter = 1;
        real_sampling_output_ant=[];
        imaginary_sampling_output_ant=[];
        sampled_h_ant = [];
       
        while(l_counter<=L)
            sampling_counter = os;
            real_sampling_output = [];
            imaginary_sampling_output = [];
            time_sampling = [];

            while(sampling_counter <= size(re_RX_output,2)-(floor(os/2)))
                real_sampling_output = [real_sampling_output re_RX_output(l_counter,sampling_counter)];
                imaginary_sampling_output = [imaginary_sampling_output im_RX_output(l_counter,sampling_counter)];
                sampling_counter = sampling_counter + os;
            end
          
            h_array = [];
            count = 1;
            index = 1;
            while(count <= size(real_sampling_output,2))
                h_array = [h_array h_ant(l_counter,index)];
                count = count + 1;
                index = index + os;
            end
            
            real_sampling_output_ant=[real_sampling_output_ant;real_sampling_output];
            imaginary_sampling_output_ant=[imaginary_sampling_output_ant;imaginary_sampling_output];
            sampled_h_ant = [sampled_h_ant;h_array];
            l_counter=l_counter+1;
        end
        
        y = real_sampling_output_ant + imaginary_sampling_output_ant * im_j;
           
        y_eq = [];
        count = 1;
        y_new = [];
        n_h_squared = [];
        while(count <= size(y,2))
            y_temp = y(:,count);
            h_temp = sampled_h_ant(:,count);
            
            mean_h = sum(h_temp)/length(h_temp);
            n_h_squared =[n_h_squared (norm(mean_h))^2];
            
            y_H = conj(transpose(h_temp));
            h_norm_squared = norm(h_temp)^2;
            y_new =[y_new ((y_H*y_temp)/h_norm_squared)];
            count = count + 1;
        end
        
        real_sampling_output = real(y_new);
        imaginary_sampling_output = imag(y_new);
    
    
    
    elseif(MIMO == 1)   %FOR THE MIMO TECHNIQUE 
        l_counter = 1;
        real_sampling_output_ant=[];
        imaginary_sampling_output_ant=[];
        sampled_h_ant = [];
       
        while(l_counter<=L)
            sampling_counter = os;
            real_sampling_output = [];
            imaginary_sampling_output = [];
            time_sampling = [];

            while(sampling_counter <= size(re_RX_output,2)-(floor(os/2)))
                real_sampling_output = [real_sampling_output re_RX_output(l_counter,sampling_counter)];
                imaginary_sampling_output = [imaginary_sampling_output im_RX_output(l_counter,sampling_counter)];
                sampling_counter = sampling_counter + os;
            end
                   
            
            real_sampling_output_ant=[real_sampling_output_ant;real_sampling_output];
            imaginary_sampling_output_ant=[imaginary_sampling_output_ant;imaginary_sampling_output];    

            l_counter=l_counter+1;
        end
        
        
        h_counter = 1;
        while(h_counter<=L*L)
            h_array = [];
            count = 1;
            index = 1;
            while(count <= size(real_sampling_output,2))
                h_array = [h_array h_ant(h_counter,index)];
                count = count + 1;
                index = index + os;
            end
                      
            sampled_h_ant = [sampled_h_ant;h_array];
            h_counter = h_counter + 1;
        end
        
        
        y = real_sampling_output_ant + imaginary_sampling_output_ant * im_j;
    
        y_eq = [];
        count = 1;
        y_new = [];
       
        while(count <= size(y,2))
            H_array = [];
            h_index_i = 1;
            while(h_index_i<=L)
                H_array_temp_j = [];
                h_index_j = 1;
                while(h_index_j<=L)
                    index = ((L*h_index_i)-L)+h_index_j;
                    H_array_temp_j = [H_array_temp_j sampled_h_ant(index,count)];                 
                    h_index_j = h_index_j + 1;
                end
                H_array = [H_array;H_array_temp_j];
                h_index_i = h_index_i + 1 ;
            end
            
            H_array_herm = conj(transpose(H_array));
            y_coef = (inv(H_array_herm*H_array))*H_array_herm;
            y_temp = y_coef*y(:,count);
           
            inter_counter = 1;
            while(inter_counter<=L)
                temp  = y_temp(inter_counter);
                y_new = [y_new temp];
                inter_counter = inter_counter+1;
            end
          
            count = count + 1;
        end
        
        real_sampling_output = real(y_new);
        imaginary_sampling_output = imag(y_new);
           
    end
  
    subplot (2, 3, 6);
    stem(time_start,real_sampling_output,"filled");
    title ("Real part of sample filter output");
    
    figure(2)
    set(figure(2), 'Position', [(1.3*s_width/2) (s_height/2) (s_width/3) (2*s_height/5)])
    scatter(real_sampling_output,imaginary_sampling_output,"filled");
    title ("Scatter plot of sampling output");
    waitforbuttonpress();
    clf(figure(1),'reset')
    clf(figure(2),'reset')
    if(spect_flag==1)
        clf(figure(4),'reset')
        clf(figure(5),'reset')
    end 
    figure(1)
    title("Computations...")
    figure(2)
    title("Computations...")
    
    
    %%%%%%%%%%%%%%%%
    %%  DECODING  %%
    %%%%%%%%%%%%%%%%
    
    output_symbols = real_sampling_output + im_j*imaginary_sampling_output;
       
    if(QAM_16 == 1)
        dataSymbolsOutG = myqamdemodulator(output_symbols,M_QAM,Eb);
        dataOutMatrixG = dec2bin(dataSymbolsOutG,k_QAM);
        dataOutG = dataOutMatrixG(:);
        final_bit_array = str2num(dataOutG);
        
    else
        output_symbols_angles = arg(output_symbols);
    
        target_angles = [];             #target values are the ideal angles of points for each class
        for i=0:m-1
            target_angles = [target_angles arg(exp(im_j*i*pi/(m/2) + im_j*(pi/m)))];
        end
        
        #in order to classify our transmitted symbols,compare the angle of each symbol with every of the ideal angles
        final_bit_array = [];
        i=1;
        while(i<=num_of_symbols)
            min_dist =1000;
            for k=0:m-1
                dist = abs(output_symbols_angles(i)-target_angles(k+1));   
                if(dist<min_dist)
                    min_dist = dist;
                    min_k = k;
                end
            end
            i=i+1;
            sub_str = dec2bin(min_k,log2(m));
            for k=1:log2(m)
                final_bit_array = [final_bit_array str2num(sub_str(1,k))];
            end
        end 
    end
    
    
    
    %%%%%%%%%%%%%%%%%%
    %%  ERROR PLOT  %%
    %%%%%%%%%%%%%%%%%%
    
    %find the number of error bits in comparison with the initial bit array
    errors = 0;
    if(QAM_16 == 0)
        Nsize =size(final_bit_array,2);
    else
        Nsize =size(final_bit_array,1);
    end
    
    
    %calculate goodput
    if(GOODPUT_CALC==1)
        correct_packets = 0;
        packet_size = Rb*Tpacket;
        
        num_of_packets = Nsize/packet_size;
        for i=1 : packet_size : Nsize
            if((i+packet_size-1)<=(Nsize))
                cur_no_of_errors = sum(abs(final_bit_array(i:(i+packet_size-1)) - input_bits(i:(i+packet_size-1))));
            else
                cur_no_of_errors = sum(abs(final_bit_array(i:(i+Nsize-i)) - input_bits(i:(i+Nsize-i))));
            end
      
            if(cur_no_of_errors == 0)
                correct_packets += 1;
            end
        end
      
        curg = correct_packets/num_of_packets;
        if(MIMO == 1)
              curg=curg*2;
        end
        goodput = [goodput curg];  
    end
    
    for i=1:Nsize
        if(xor(input_bits(i),final_bit_array(i))==1)
            errors = errors +1;     
        end
    end
    
    errors_array = [errors_array errors];
    cur_BER = errors/Nsize;
    
    BER = [BER cur_BER];
    
    SNRSlin = [SNRSlin SNRlin];
    if(MIMO==0 && (OFDM==0))
        mn = [mn sum(n_h_squared)/length(n_h_squared)];
    end
    snr_counter = snr_counter+1;
end



%For the theoretical BER
if(MIMO == 1 && OFDM==0)
    count = 1;
    size(SNRSlin,2)
    while(count <= size(SNRSlin,2))
        gamma = SNRSlin(count);
        BERavg = [BERavg (0.5*(1-sqrt((gamma)/(1+gamma))))];
        count = count + 1;
    end
elseif(OFDM==0)
    count = 1;
    final_mean = sum(mn)/length(mn);
    while(count <= length(mn))
        gamma = SNRSlin(count)*final_mean;
        BERavg = [BERavg (0.5*(1-sqrt((gamma)/(1+gamma))))];
        count = count + 1;
    end
end
    
BER
if(GOODPUT_CALC == 1)
    goodput
end
%Create the BER to SNR plot

if(OFDM == 1)
    f3 = figure(3);
    set(f3, 'Position', [(1.5*s_width/5) 0 (2*s_width/5) (1.6*s_height/5)])
    t=SNRdB(1:end);
    d=log(theoretical_BER(1:end));
    y1 = semilogy(SNRdB,theoretical_BER,"-og",'linewidth',2);
    hold on;
    y2 = semilogy(SNRdB,BER,"-or",'linewidth',2);
    title("BER vs SNR(dB)");
    ylabel("BER");
    xlabel("SNR(dB)");
    hh = legend ({"Theoretical AWGN channel BER"}, " BER with OFDM");
    legend(hh, "location", "northoutside")

elseif((d_flag==1) && (MIMO == 0))
    f3 = figure(3);
    set(f3, 'Position', [(1.5*s_width/5) 0 (2*s_width/5) (1.6*s_height/5)])

    t=SNRdB(1:end);
    d=log(theoretical_BER(1:end));
    y1 = semilogy(SNRdB,theoretical_BER,"-og",'linewidth',2);
    hold on;
    y2 = semilogy(SNRdB,BER,"-or",'linewidth',2);
    ylabel("BER");
    xlabel("SNR(dB)");
    hh = legend ({"Theoretical AWGN channel BER"}, " BER with equalization/channel inversion using MRC method");
    legend(hh, "location", "northoutside")

elseif((d_flag == 0)  && (MIMO == 0))
    f3 = figure(3);
    set(f3, 'Position', [(1.5*s_width/5) 0 (2*s_width/5) (1.6*s_height/5)])

    t=SNRdB(1:end);
    d=log(theoretical_BER(1:end));
    y1 = semilogy(SNRdB,theoretical_BER,"-og",'linewidth',2);
    hold on;
    y2 = semilogy(SNRdB,BER,"-or",'linewidth',2);
    hold on;
    y3 = semilogy(SNRdB,BERavg,"-ob",'linewidth',2);
    ylabel("BER");
    xlabel("SNR(dB)");
    hh = legend ({"Theoretical AWGN channel BER"}, " BER with equalization/channel inversion","Average theoretical BER with equalization/channel inversion");
    legend(hh, "location", "northoutside")

elseif(MIMO == 1)

    f3 = figure(3);
    set(f3, 'Position', [(1.5*s_width/5) 0 (2*s_width/5) (1.6*s_height/5)])

    t=SNRdB(1:end);
    d=log(theoretical_BER(1:end));
    y1 = semilogy(SNRdB,theoretical_BER,"-og",'linewidth',2);
    hold on;
    y2 = semilogy(SNRdB,BER,"-or",'linewidth',2);
    hold on;
    y3 = semilogy(SNRdB,BERavg,"-ob",'linewidth',2);
    ylabel("BER");
    xlabel("SNR(dB)");
    hh = legend ({"Theoretical AWGN channel BER"}, " BER with MIMO","Average theoretical BER of MIMO");
    legend(hh, "location", "northoutside")

end

waitforbuttonpress();

close(figure(1));
close(figure(2));
close(figure(3));   

if(spect_flag==1)
    close(figure(4));
    close(figure(5));
end
     