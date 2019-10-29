%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Task 1         %%
close all; clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup (part A)
beta = 4;   % Number of bits per symbol--BPSK,QPSK ...
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

xencoded =[g0'; g1'];       % Encoded x
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
SNRdB = 50;       % SNR
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


%% Multiple SNR 

% Xup = mapping(xinterleaved,Mapping);
% SNRdB = [-20:1:20];
% % Adding noise
% var = Eup/(2*10^(SNRdB/10));    % var^2 calculated
% noise = randn(1,Nn)*sqrt(var);  % AWGN
% Signal = Xup + noise; 
% % demapping
% 
% % de-Interweaving
% Xallinter = zeros(length(SNRdB),length(Xup));
% for u = 1:length(SNRdB)
%     Xd = Xalld(i,1:length(Xup));
%     for i = 1:N:length(Xd)
%         tmp = Xd(i:i+N-1);
%         tmp = reshape(tmp,16,16).';
%         tmp = reshape(tmp,1,N);
%         I(i:i+N-1) = tmp;
%     end
%     Xallinter(i,1:length(Xup)) = I;
% end
% 
% vitdec
% decodeds = zeros(length(SNRdB),length(decoded));
% for i = 1:length(SNRdB)
%     tr = poly2trellis([m + 1], [g(1) g(2)]);
%     tracebacklength = 5*m;
%     decodeds(i,1:length(decoded)) = vitdec(Xallinter(i,1:length(Xup)), tr, tracebacklength, 'term', 'hard'); 
% end








