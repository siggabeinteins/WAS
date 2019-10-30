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

xencoded =[g0; g1];       % Encoded x
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
iterqpsk = 300;
iterbpsk = 500;
iter16 = 200;
iter64 = 200;
    %% QPSK Errsqpsk
SNRdB = [-20:1:20];
Errsqpsk = zeros(1,length(SNRdB));
fprintf('QPSK ');
for t=1:length(SNRdB)
    errors = zeros(1,length(SNRdB));
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
    errors = zeros(1,length(SNRdB));
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
    errors = zeros(1,length(SNRdB));
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
    errors = zeros(1,length(SNRdB));
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
fprintf('\n');
    %% Plotting for part K
figure
semilogy(SNRdB,Errsbpsk)
hold on
semilogy(SNRdB,Errsqpsk,'r')
semilogy(SNRdB,Errs16,'g')
semilogy(SNRdB,Errs64,'c')
legend('BPSK','QPSK','16-QAM','64-QAM')

%% Uncoded transmission (Part L)
    %% QPSK Errsqpsk2
Errsqpsk2 = zeros(1,length(SNRdB));
fprintf('\nQPSK-uncoded ');
beta = 2;
Nb = N;
Mapping = 'QPSK';
for t=1:length(SNRdB)
    errors = zeros(1,length(SNRdB));
    totalerror = 0;
    SNRdBi = SNRdB(t);
    for u = 1:iterqpsk
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
    Errsqpsk2(t) = totalerror/iterqpsk;
    fprintf('%d ',t);
end
    %% Plotting for part L
figure
semilogy(SNRdB,Errsqpsk)
hold on
semilogy(SNRdB,Errsqpsk2)
legend('Encoded','Non-encoded')









