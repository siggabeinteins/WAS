%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Task 1         %%
close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup (part A)
beta = 2;   % Number of bits per symbol--BPSK,QPSK ...
N = 256;    % Length of one block of symbols
m = 3;      % Time delay, dependant of number of conv
R = 1/2;    % Rate of code
Nb = N*beta*R-m; % Length of input bit sequence

% Generation of random bits
x = randi([0,1],1,Nb);

%% Encoding (part B)
switch m
    case 2
        co0 = [1 0 1];
        co1 = [1 1 1];
        g = [5 7];
    case 3
        co0 = [1 0 1 1];
        co1 = [1 1 1 1];
        g = [13 17];
    case 4
        co0 = [1 0 1 1 1];
        co1 = [1 1 0 0 1];
        g = [27 31];
    case 5
        co0 = [1 0 1 0 1 1];
        co1 = [1 1 1 1 0 1];
        g = [53 75];
    case 6
        co0 = [1 0 0 1 1 1 1];
        co1 = [1 1 0 1 1 0 1];
        g = [117 155];
end

g0 = conv(x,co0);           % Convolving with the convolution matrix
g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
g1 = conv(x,co1);           % Convolving with the convolution matrix
g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix

xencoded =[g0; g1];         % Encoded x
xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
%xencoded = 1:1536;         % Good for visualization

%% Drawing diff m (Part C)
% here I drew pictures and will just displaying them

%% Interleaving (Part D)
Blocklength = N;        % Size of block interleaver
I = zeros(1,length(xencoded));    % Initalize the Interleaved signal

% The for loop splits the Xencoded signal into blocks, interleaves the
% block and then appends it to the interleaved output signal
for i = 1:N:length(xencoded)
    tmp = xencoded(i:i+N-1);
    tmp = reshape(tmp.',16,16).';
    tmp = reshape(tmp,1,N);
    I(i:i+N-1) = tmp;
end
% Here we get a signal that where the blocks have been interleaved
xinterleaved = I;

%% Mapping (Part E)
%----------------------
Mapping = 'QPSK';       % Change to choose mapping
%----------------------
Xup = mapping(xinterleaved,Mapping);

%% AWGN added (Part F)
%----------------------
SNRdB = 10;       % SNR
%----------------------
Eup = 0; Nn = length(Xup);      % Initialization
% Calculation of total power
for k = 1:Nn
    Eup = Eup+abs(Xup(k))^2;
end
Eup = Eup/Nn;                   % Eup
var = Eup/(2*10^(SNRdB/10));    % var^2 calculated
noise = randn(1,Nn)*sqrt(var);  % AWGN
Signal = Xup + noise;           % Signal with noise

%% Demapping (Part G)
Xdown = demapping(Signal,Mapping);  % Signal Demapped 

%% De-interleaving(Part H)
I = zeros(1,length(Xdown));         % Initialization
for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
    tmp = Xdown(i:i+N-1);
    tmp = reshape(tmp,16,16).';
    tmp = reshape(tmp,1,N);
    I(i:i+N-1) = tmp;
end
Xdinterleaved = I;

%% Decoding Viterbi(Part I)
% Generate the trellis from the polynomials
tr = poly2trellis([m + 1], [g(1) g(2)]);
% Set the traceback length. If the code rate is 1/2, a typical value for tracebacklength is
% about five times m.
tracebacklength = 5*m;
% Call the vitdec function.
decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');   % coded is the coded input
                                                                % bits.
% The output of vitdec is the decoded block of bits including the trailing bits.


%% BER for this run(Part J)
Nerr = 0;                           % Initialization
for i=1:Nb                          % Tally up all errors 
    if decoded(i)~=x(i)
        Nerr=Nerr+1;
    end
end
Errorrate = 100*Nerr/Nb;            % Error rate calculated


%% Multiple SNR (Part K)
% Here we can decide the iterations for each mapping type
iterqpsk = 100;
iterbpsk = 200;
iter16 = 100;
iter64 = 100;
SNRdB = -20:1:20;
    %% QPSK Errsqpsk
Errsqpsk = zeros(1,length(SNRdB));
fprintf('QPSK ');
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iterqpsk
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsqpsk(t) = totalerror/iterqpsk;
    fprintf('%d ',t);
end
    %% BPSK Errsbpsk
Errsbpsk = zeros(1,length(SNRdB));
fprintf('\nBPSK ');
beta = 1;
Nb = N*beta*R-m;
Mapping = 'BPSK';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iterbpsk
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsbpsk(t) = totalerror/iterbpsk;
    fprintf('%d ',t);
end
    %% 16-QAM Errs16
Errs16 = zeros(1,length(SNRdB));
fprintf('\n16-QAM ');
beta = 4;
Nb = N*beta*R-m;
Mapping = '16-QAM';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iter16
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs16(t) = totalerror/iter16;
    fprintf('%d ',t);
end
    %% 64-QAM Errs64
Errs64 = zeros(1,length(SNRdB));
fprintf('\n64-QAM ');
beta = 6;
Nb = N*beta*R-m;
Mapping = '64-QAM';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iter64
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs64(t) = totalerror/iter64;
    fprintf('%d ',t);
end
    %% Plotting for part K
figure
semilogy(SNRdB,Errsbpsk)
hold on
semilogy(SNRdB,Errsqpsk,'r')
semilogy(SNRdB,Errs16,'g')
semilogy(SNRdB,Errs64,'c')
legend('BPSK','QPSK','16-QAM','64-QAM')
ylim([0.2 100])
xlabel('SNR [dB]');
ylabel('Errorrate[%]');
fprintf('\n');

%% Uncoded transmission (Part L)
iterqpsk2 = 100;
iterbpsk2 = 200;
iter162 = 100;
iter642 = 100;
    %% QPSK Errsqpsk2
Errsqpsk2 = zeros(1,length(SNRdB));
fprintf('\nQPSK-uncoded ');
beta = 2;
Nb = N*beta;
Mapping = 'QPSK';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iterqpsk2
        x = randi([0,1],1,Nb);
%         g0 = conv(x,co0);           % Convolving with the convolution matrix
%         g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
%         g1 = conv(x,co1);           % Convolving with the convolution matrix
%         g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
%         xencoded =[g0; g1];         % Encoded x
%         xencoded = x;    % Interleaving g0 and g1 to make the encoded x
%         I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
%         for i = 1:N:length(xencoded)
%             tmp = xencoded(i:i+N-1);
%             tmp = reshape(tmp.',16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         xinterleaved = I;
        Xup = mapping(x,Mapping);
        Nn = length(Xup);
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
%         I = zeros(1,length(Xdown));         % Initialization
%         for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
%             tmp = Xdown(i:i+N-1);
%             tmp = reshape(tmp,16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         Xdinterleaved = I;
%         tr = poly2trellis([m + 1], [g(1) g(2)]);
%         tracebacklength = 5*m;
%         decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        decoded = Xdown;
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsqpsk2(t) = totalerror/iterqpsk2;
    fprintf('%d ',t);
end
    %% BPSK Errsbpsk2
Errsbpsk2 = zeros(1,length(SNRdB));
fprintf('\nBPSK-uncoded ');
beta = 1;
Nb = N*beta;
Mapping = 'BPSK';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iterbpsk
        x = randi([0,1],1,Nb);
%         g0 = conv(x,co0);           % Convolving with the convolution matrix
%         g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
%         g1 = conv(x,co1);           % Convolving with the convolution matrix
%         g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
%         xencoded =[g0; g1];         % Encoded x
%         xencoded = x;    % Interleaving g0 and g1 to make the encoded x
%         I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
%         for i = 1:N:length(xencoded)
%             tmp = xencoded(i:i+N-1);
%             tmp = reshape(tmp.',16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         xinterleaved = I;
        Xup = mapping(x,Mapping);
        Nn = length(Xup);
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
%         I = zeros(1,length(Xdown));         % Initialization
%         for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
%             tmp = Xdown(i:i+N-1);
%             tmp = reshape(tmp,16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         Xdinterleaved = I;
%         tr = poly2trellis([m + 1], [g(1) g(2)]);
%         tracebacklength = 5*m;
%         decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        decoded = Xdown;
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsbpsk2(t) = totalerror/iterbpsk;
    fprintf('%d ',t);
end
    %% 16-QM Errs162
Errs162 = zeros(1,length(SNRdB));
fprintf('\n16-QAM-uncoded ');
beta = 4;
Nb = N*beta;
Mapping = '16-QAM';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iter162
        x = randi([0,1],1,Nb);
%         g0 = conv(x,co0);           % Convolving with the convolution matrix
%         g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
%         g1 = conv(x,co1);           % Convolving with the convolution matrix
%         g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
%         xencoded =[g0; g1];         % Encoded x
%         xencoded = x;    % Interleaving g0 and g1 to make the encoded x
%         I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
%         for i = 1:N:length(xencoded)
%             tmp = xencoded(i:i+N-1);
%             tmp = reshape(tmp.',16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         xinterleaved = I;
        Xup = mapping(x,Mapping);
        Nn = length(Xup);
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
%         I = zeros(1,length(Xdown));         % Initialization
%         for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
%             tmp = Xdown(i:i+N-1);
%             tmp = reshape(tmp,16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         Xdinterleaved = I;
%         tr = poly2trellis([m + 1], [g(1) g(2)]);
%         tracebacklength = 5*m;
%         decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        decoded = Xdown;
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs162(t) = totalerror/iter162;
    fprintf('%d ',t);
end
    %% 64-QM Errs642
Errs642 = zeros(1,length(SNRdB));
fprintf('\n64-QAM-uncoded ');
beta = 6;
Nb = N*beta;
Mapping = '64-QAM';
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iter642
        x = randi([0,1],1,Nb);
%         g0 = conv(x,co0);           % Convolving with the convolution matrix
%         g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
%         g1 = conv(x,co1);           % Convolving with the convolution matrix
%         g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
%         xencoded =[g0; g1];         % Encoded x
%         xencoded = x;    % Interleaving g0 and g1 to make the encoded x
%         I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
%         for i = 1:N:length(xencoded)
%             tmp = xencoded(i:i+N-1);
%             tmp = reshape(tmp.',16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         xinterleaved = I;
        Xup = mapping(x,Mapping);
        Nn = length(Xup);
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
%         I = zeros(1,length(Xdown));         % Initialization
%         for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
%             tmp = Xdown(i:i+N-1);
%             tmp = reshape(tmp,16,16).';
%             tmp = reshape(tmp,1,N);
%             I(i:i+N-1) = tmp;
%         end
%         Xdinterleaved = I;
%         tr = poly2trellis([m + 1], [g(1) g(2)]);
%         tracebacklength = 5*m;
%         decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        decoded = Xdown;
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs642(t) = totalerror/iter642;
    fprintf('%d ',t);
end
    %% Plotting for part L
figure
semilogy(SNRdB,Errsqpsk)
hold on
semilogy(SNRdB,Errsqpsk2)
legend('Encoded','Non-encoded')
title('QPSK')
xlabel('SNR [dB]');
ylabel('Errorrate[%]');

figure
semilogy(SNRdB,Errsbpsk)
hold on
semilogy(SNRdB,Errsbpsk2)
legend('Encoded','Non-encoded')
title('BPSK')
xlabel('SNR [dB]');
ylabel('Errorrate[%]');

figure
semilogy(SNRdB,Errs16)
hold on
semilogy(SNRdB,Errs162)
legend('Encoded','Non-encoded')
title('16-QAM')
xlabel('SNR [dB]');
ylabel('Errorrate[%]');

figure
semilogy(SNRdB,Errs64)
hold on
semilogy(SNRdB,Errs642)
legend('Encoded','Non-encoded')
title('64-QAM')
xlabel('SNR [dB]');
ylabel('Errorrate[%]');
fprintf('\n');

%% Different m´s (Part M)
biter = 200;
ite64 = 100;
    %% Begin with m = 2
m = 2;
fprintf('\nm=2 ');
co0 = [1 0 1];
co1 = [1 1 1];
g = [5 7];
        %% First calculate BPSK
Errsbpskm2 = zeros(1,length(SNRdB));
fprintf('\nBPSK 2');
Mapping = 'BPSK';
beta = 1;   % Number of bits per symbol--BPSK,QPSK ...
Nb = N*beta*R-m; % Length of input bit sequence
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:biter
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsbpskm2(t) = totalerror/biter;
    fprintf('%d ',t);
end
        %% Then calculate 64-QAM 
Errs64m2 = zeros(1,length(SNRdB));
fprintf('\n64-QAM 2');
Mapping = '64-QAM';
beta = 6;   % Number of bits per symbol--BPSK,QPSK ...
Nb = N*beta*R-m; % Length of input bit sequence
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:ite64
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs64m2(t) = totalerror/ite64;
    fprintf('%d ',t);
end        
    %% m = 4
m = 4;
fprintf('\nm=4 ');
co0 = [1 0 1 1 1];
co1 = [1 1 0 0 1];
g = [27 31];
        %% First calculate BPSK
Errsbpskm4 = zeros(1,length(SNRdB));
fprintf('\nBPSK 4');
Mapping = 'BPSK';
beta = 1;   % Number of bits per symbol--BPSK,QPSK ...
Nb = N*beta*R-m; % Length of input bit sequence
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:biter
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsbpskm4(t) = totalerror/biter;
    fprintf('%d ',t);
end
        %% Then calculate 64-QAM 
Errs64m4 = zeros(1,length(SNRdB));
fprintf('\n64-QAM 4');
Mapping = '64-QAM';
beta = 6;   % Number of bits per symbol--BPSK,QPSK ...
Nb = N*beta*R-m; % Length of input bit sequence
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:ite64
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs64m4(t) = totalerror/ite64;
    fprintf('%d ',t);
end        
    %% m = 6
m = 6;
fprintf('\nm=6 ');
co0 = [1 0 0 1 1 1 1];
co1 = [1 1 0 1 1 0 1];
g = [117 155];
        %% First calculate BPSK
Errsbpskm6 = zeros(1,length(SNRdB));
fprintf('\nBPSK 6');
Mapping = 'BPSK';
beta = 1;   % Number of bits per symbol--BPSK,QPSK ...
Nb = N*beta*R-m; % Length of input bit sequence
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:biter
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errsbpskm6(t) = totalerror/biter;
    fprintf('%d ',t);
end
        %% Then calculate 64-QAM 
Errs64m6 = zeros(1,length(SNRdB));
fprintf('\n64-QAM 6');
Mapping = '64-QAM';
beta = 6;   % Number of bits per symbol--BPSK,QPSK ...
Nb = N*beta*R-m; % Length of input bit sequence
for t=1:length(SNRdB)
    
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:ite64
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Xup = mapping(xinterleaved,Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errs64m6(t) = totalerror/ite64;
    fprintf('%d ',t);
end        
        %% Plotting for part L
figure
semilogy(SNRdB,Errsbpskm2)
hold on
semilogy(SNRdB,Errsbpskm4,'r')
semilogy(SNRdB,Errsbpskm6,'g')
semilogy(SNRdB,Errs64m2,'-.b')
semilogy(SNRdB,Errs64m4,'-.r')
semilogy(SNRdB,Errs64m6,'-.g')
title('Difference from m = 2,4,6')
legend('BPSK,m=2','BPSK,m=4','BPSK,m=6','64-QAM,m=2','64-QAM,m=4','64-QAM,m=6')
xlabel('SNR [dB]');
ylabel('Errorrate[%]');

%% Trying out (Part K)
fprintf('\n');
partkiters = 500;
co0 = [1 0 1 1 1];
co1 = [1 1 0 0 1];
g = [27 31];
beta = 4;
m = 4;
Mapping = '16-QAM';
SNRdB = -10:2:30;
Nb = N*beta*R-m; % Length of input bit sequence
    %% Burst = 1
b = 1;
Errspartk1 = zeros(1,length(SNRdB));
fprintf('\n16-QAM part b=1 ');
for t=1:length(SNRdB) 
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:partkiters
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Z = ones(1,1024);
        for i = 1:256:1024
            li = randi([1 255],1,1);
            Z(i+li-b+1:i+li) = zeros(1,b);
        end
        Xup = mapping(xinterleaved.*Z(1:length(xinterleaved)),Mapping);    
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errspartk1(t) = totalerror/partkiters;
    fprintf('%d ',t);
end 
    %% Burst = 2
b = 2;
Errspartk2 = zeros(1,length(SNRdB));
fprintf('\n16-QAM part b=2 ');
for t=1:length(SNRdB) 
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:partkiters
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Z = ones(1,1024);
        for i = 1:256:1024
            li = randi([1 255],1,1);
            Z(i+li-b+1:i+li) = zeros(1,b);
        end
        Xup = mapping(xinterleaved.*Z(1:length(xinterleaved)),Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errspartk2(t) = totalerror/partkiters;
    fprintf('%d ',t);
end 
    %% Burst = 3
b = 3;
Errspartk3 = zeros(1,length(SNRdB));
fprintf('\n16-QAM part b=3 ');
for t=1:length(SNRdB) 
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:partkiters
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Z = ones(1,1024);
        for i = 1:256:1024
            li = randi([b+1 255],1,1);
            Z(i+li-b+1:i+li) = zeros(1,b);
        end
        Xup = mapping(xinterleaved.*Z(1:length(xinterleaved)),Mapping);    
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errspartk3(t) = totalerror/partkiters;
    fprintf('%d ',t);
end 
    %% Burst = 4
b = 4;
Errspartk4 = zeros(1,length(SNRdB));
fprintf('\n16-QAM part b=4 ');
for t=1:length(SNRdB) 
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:partkiters
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Z = ones(1,1024);
        for i = 1:256:1024
            li = randi([b+1 255],1,1);
            Z(i+li-b+1:i+li) = zeros(1,b);
        end
        Xup = mapping(xinterleaved.*Z(1:length(xinterleaved)),Mapping);    
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errspartk4(t) = totalerror/partkiters;
    fprintf('%d ',t);
end 
    %% Burst = 5
b = 5;
Errspartk5 = zeros(1,length(SNRdB));
fprintf('\n16-QAM part b=5 ');
for t=1:length(SNRdB) 
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:partkiters
        x = randi([0,1],1,Nb);
        g0 = conv(x,co0);           % Convolving with the convolution matrix
        g0 = mod(g0,2);             % Taking modulo 2 of the whole matrix
        g1 = conv(x,co1);           % Convolving with the convolution matrix
        g1 = mod(g1,2);             % Taking modulo 2 of the whole matrix
        xencoded =[g0; g1];         % Encoded x
        xencoded = xencoded(:)';    % Interleaving g0 and g1 to make the encoded x
        I = zeros(1,length(xencoded));    % Initalize the Interleaved signal
        for i = 1:N:length(xencoded)
            tmp = xencoded(i:i+N-1);
            tmp = reshape(tmp.',16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        xinterleaved = I;
        Z = ones(1,1024);
        for i = 1:256:1024
            li = randi([b+1 255],1,1);
            Z(i+li-b+1:i+li) = zeros(1,b);
        end
        Xup = mapping(xinterleaved.*Z(1:length(xinterleaved)),Mapping);   
        var = Eup/(2*10^(SNRdBi/10));    % var^2 calculated
        noise = randn(1,Nn)*sqrt(var);  % AWGN
        Signal = Xup + noise; 
        Xdown = demapping(Signal,Mapping);
        I = zeros(1,length(Xdown));         % Initialization
        for i = 1:N:length(Xdown)           % The incomming signal is deinterleaved
            tmp = Xdown(i:i+N-1);
            tmp = reshape(tmp,16,16).';
            tmp = reshape(tmp,1,N);
            I(i:i+N-1) = tmp;
        end
        Xdinterleaved = I;
        tr = poly2trellis([m + 1], [g(1) g(2)]);
        tracebacklength = 5*m;
        decoded = vitdec(Xdinterleaved, tr, tracebacklength, 'term', 'hard');
        Nerr = 0;                           % Initialization
        for i=1:Nb                          % Tally up all errors 
            if decoded(i)~=x(i)
                Nerr=Nerr+1;
            end
        end
        Errorrate = 100*Nerr/Nb;            % Error rate calculated
        totalerror = totalerror + Errorrate;
    end
    Errspartk5(t) = totalerror/partkiters;
    fprintf('%d ',t);
end 
    %% Plotting for part K
figure
semilogy(SNRdB,Errspartk1)
hold on
semilogy(SNRdB,Errspartk2,'r')
semilogy(SNRdB,Errspartk3,'g')
semilogy(SNRdB,Errspartk4,'c')
semilogy(SNRdB,Errspartk5,'y')
title('16-QAM, m = 4 , different burst rates of error')
legend('Burst = 1','Burst = 2','Burst = 3','Burst = 4','Burst = 5')
xlabel('SNR [dB]');
ylabel('Errorrate[%]');