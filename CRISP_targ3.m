%code to 1) extract endogenous candidate Cas9 binding sites from a target sequence, and 2) screen them for orthogonality to context sequence (ex.: genome) 

%***load sequence context sequence (ex. genome), and target sequence from
%local file
[HeadC, SeqC] = fastaread('E.coli-MG1655-genome.fasta');
%[HeadS, SeqT] = fastaread(File);

%***alternatively load directly form Genbank
%S = getgenbank('M10051');


%generate (reverse) complement
SeqC_C = seqrcomplement(SeqC);
SeqC_RC = seqrcomplement(SeqC);

%convert sequence to vector of integers where A=1, G=3, C=2, T=4
SeqC_v = nt2int(SeqC);
SeqC_C_v = nt2int(SeqC_C);
SeqC_RC_v = nt2int(SeqC_RC);


%***alternative to importing, isolate target sequence from context sequence
%%eventually SeqT_v = getgenbank('gene'); ...

%**define start and length of target sequence
start = 4.637613e6; len = 716;
%in this case, we are targeting the ArcA gene

%extract target sequence (and its reverse complement) from context sequence
SeqT_v = SeqC_v(start:start+len-1); 
SeqT_C_v = SeqC_C_v(start:start+len-1); 
SeqT_RC_v = fliplr(SeqT_C_v);


%**length of gRNA subsequence 
lseq=23;

%***scan target sequences to find GG positions and output a library of their
%associated gRNA candidates 
[ sSeqT_v, sSeqT_RC_v ] = GG_library2( SeqT_v, SeqT_RC_v, lseq); 

%glue and display proportion subsequence candidates
sSeqT_glued = [sSeqT_v ; sSeqT_RC_v];
lf = length(sSeqT_v);
disp(['number of target targeting candidates in forward library = ' int2str(lf)])
disp(['number of target targeting candidates in reverse complement library = ' int2str(length(sSeqT_RC_v))])



%cut out target sequence (minus 22 nt (lseq-1) overhangs) from context sequence
SeqCmT_v_LHS = SeqC_v(1:start+lseq-1); 
SeqCmT_v_RHS = SeqC_v(start+len-(lseq-1):end); 
SeqCmT_RC_v_LHS = fliplr(SeqC_C_v(1:start+lseq-1)); 
SeqCmT_RC_v_RHS = fliplr(SeqC_C_v(start+len-(lseq-1):end)); 
SeqCmT_glued_LHS = [SeqCmT_v_LHS ; SeqCmT_RC_v_LHS];
SeqCmT_glued_RHS = [SeqCmT_v_RHS ; SeqCmT_RC_v_RHS];


%***with the library of small targeting sequences, scan through context sequence and calculate weighted hamming distance at
%each position (with mostly-excised target sequence) 
[ m_L, m_R ]=scan_SeqC3( sSeqT_glued, lf, SeqCmT_glued_LHS, SeqCmT_glued_RHS, lseq ); 

%concatenate scores for left and right hand sides
m1 = [m_L m_R];

%find highest score from all translations across the context sequence
% make that score representative of subsequences ability to bind somewhere
% on the context sequence
maxscores = max(m1');

scatter(1:length(sSeqT_glued(:,lseq)), maxscores)

% identify subsequence that is most orthogonal to any part of the context sequence. 
[m2, m2i] = min(maxscores);

%output result
disp(['best (most orthogonal to context) target subsequence is #' int2str(m2i) ' in sSeqT_glued'])

int2nt(sSeqT_glued(m2i,:))
