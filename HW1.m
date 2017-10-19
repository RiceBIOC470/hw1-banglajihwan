GB Comments

Prob1: 100%
Prob2:
P1:100
P2: 100
P3:100
P4:90 Could use some labels.  
P5:100 
Prob3
P1: 100
P2:100
P3:100 
Overall: 99


% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
% x = 3; y = '5'; %mixed

%your code goes here
if isnumeric(x) == 1;
    if isnumeric(y) == 1 
           z= x+y;
           disp(z) 
    
    else ischar(y) == 1; 
        z = x+str2num(y);
        disp(z) 
    end 
else
    if isnumeric(y) == 1;
           z= str2num(x)+y;
           disp(z) 
    
    else ischar(y) == 1; 
            z = str2num(x)+str2num(y);
            disp(z)  
        end 
    end
     
%output your answer
% 8 for all the cases above. 
%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length
nucleotides = ['A','T','G','C'] ;
random_num = randi(4,1,N);
rand_seq= nucleotides(random_num) 

% multiple of three! 
%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.
start = strfind(rand_seq, 'ATG');
stop_1 = strfind(rand_seq, 'TAA');
stop_2 =strfind(rand_seq, 'TGA');
stop_3 = strfind(rand_seq, 'TAG') ;
stop_all= [stop_1 stop_2 stop_3];
length=[];
seq=[];
stop_all2=sort(stop_all, 'descend');%sort from biggest to smallest
for i = stop_all2 
    for j = 1:size(start,1)
        length = i-start(j);
        %start(j)
        if rem(length,3)==0;
            seq=rand_seq(start(j):i+2);
            varA = 1; 
            break %the loop breaks once we find the difference between start and stop position divisible by 3. 
       
        else 
            length=[];
            disp('none');
            
        end 
    end 
    if varA == 1
        break 
    end 
end 
disp(seq)
%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.
con = 1; 
count = 0;
while con < 1000
    disp(con)
    N = 500; % define sequence length
    nucleotides = ['A','T','G','C'] ;
    random_num = randi(4,1,N);
    rand_seq= nucleotides(random_num); 
    start = strfind(rand_seq, 'ATG')';
    stop_1 = strfind(rand_seq, 'TAA');
    stop_2 =strfind(rand_seq, 'TGA');
    stop_3 = strfind(rand_seq, 'TAG') ;
    stop_all= [stop_1 stop_2 stop_3];
length=[];
seq=[];
stop_all2=sort(stop_all, 'descend');%sort from biggest to smallest
for i = stop_all2 
    for j = 1:size(start,1)
        length = i-start(j);  
        if rem(length,3)==0;
            seq=rand_seq(start(j):i+2);
            varA = 1; 
            break %the loop breaks once we find the difference between start and stop position divisible by 3. 
        else 
            length=[];        
        end 
    end 
    if varA == 1
        break 
    end 
end 
    seq_length=size(seq,2); 
    if seq_length <=50
        count=count+1;
    else 
    end 
    con=con+1;
end 
disp(count)
probab=1-(count/1000)
% I got 91.0% 

%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 

N=1; 
M=[ ];
P = [ ];
while N < 550
con = 1; 
count = 0;
disp(N) 
while con < 1000
    disp(con)
    % define sequence length
    nucleotides = ['A','T','G','C'] ;
    random_num = randi(4,1,N);
    rand_seq= nucleotides(random_num); 
    start = strfind(rand_seq, 'ATG')';
    stop_1 = strfind(rand_seq, 'TAA');
    stop_2 =strfind(rand_seq, 'TGA');
    stop_3 = strfind(rand_seq, 'TAG') ;
    stop_all= [stop_1 stop_2 stop_3];
length=[];
seq=[];
stop_all2=sort(stop_all, 'descend');%sort from biggest to smallest
for i = stop_all2 
    for j = 1:size(start,1)
        length = i-start(j);  
        if rem(length,3)==0;
            seq=rand_seq(start(j):i+2);
            varA = 1; 
            break %the loop breaks once we find the difference between start and stop position divisible by 3. 
        else 
            length=[];        
        end 
    end 
    if varA == 1
        break 
    end 
end 
    seq_length=size(seq,2); 
    if seq_length <=50
        count=count+1;
    else 
    end 
    con=con+1;
end 
disp(count)
probab=1-(count/1000)

M = [M,N];
P = [P,probab];
N=N+1;

end 
figure
plot(M,P)
%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this.

%When N is small we will see smaller chance of finding ORF greater than
%50bp. As N increases, the probability of finding an ORF greater than 50bp go towards 1. The curve in-between
%low and high N forms a hyperbolic curve. 


%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 
fid = fopen('qPCRdata.txt', 'r');
C = textscan(fid,'%s %s %s %s %s %f32','headerlines',2);
fclose(fid);
X=C{6}
X(isnan(X))=[]
X_1=X(1:72)%X_1 has cp data without G and H


% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 
X_2= reshape(X_1,12,6);
X_96 = [];
for i = 1:6
   X_96 =[X_96; X_2(:,i)']; %X_96 gives you the format of the well 
end 
% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

%Normalization
CpN0 = sum(X_96(1,10:12))/3;
CpNX_minus1 = [];
for i = 2:6
    CpNX = sum(X_96(i, 10:12))/3;
    CpNX_minus = CpN0 - CpNX;
    CpNX_minus1 = [CpNX_minus1; CpNX_minus]
end 

%Gene1
Cp10 = sum(X_96(1,1:3))/3;
Cp1_minus1 = [];
for i = 2:6
    Cp1 = sum(X_96(i, 1:3))/3;
    Cp1_minus = Cp10 - Cp1;
    Cp1_minus1 = [Cp1_minus1; Cp1_minus];
end 
power1 = Cp1_minus1 - CpNX_minus1;
gene1 = 2.^power1;

%Gene2
Cp20 = sum(X_96(1,4:6))/3;
Cp2_minus1 = [];
for i = 2:6
    Cp2 = sum(X_96(i, 4:6))/3;
    Cp2_minus = Cp20 - Cp2;
    Cp2_minus1 = [Cp2_minus1; Cp2_minus];
end 
power2 = Cp2_minus1 - CpNX_minus1;
gene2 = 2.^power2;

%Gene3
Cp30 = sum(X_96(1,7:9))/3;
Cp3_minus1 = [];
for i = 2:6
    Cp3 = sum(X_96(i, 7:9))/3;
    Cp3_minus = Cp30 - Cp3;
    Cp3_minus1 = [Cp3_minus1; Cp3_minus];
end 
power3 = Cp3_minus1 - CpNX_minus1;
gene3 = 2.^power3;

figure 
ax1 = subplot(2,2,1);
bar(ax1, gene1);
condition = {'conditon1', 'conditon2', 'conditon3', 'conditon4','conditon5','conditon6'}; 
set(gca,'xticklabel',condition)
ylabel ('fold-change');
title('gene1')

ax2 = subplot(2,2,2);
bar(ax2, gene2);
condition = {'conditon1', 'conditon2', 'conditon3', 'conditon4','conditon5','conditon6'} ;
set(gca,'xticklabel',condition)
title('gene2')

ax3 = subplot(2,2,3);
ylabel ('fold-change');
bar(ax3, gene3)
condition = {'conditon1', 'conditon2', 'conditon3', 'conditon4','conditon5','conditon6'} ;
set(gca,'xticklabel',condition)
ylabel ('fold-change');
title('gene3')
%row_names = {'condition1', 'condition2', 'condition3', 'condition4', 'condition5'};   
%column_names = {'gene1', 'gene2', 'gene3'}
%T = table(gene1, gene2, gene3, 'RowNames', row_names)

%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


