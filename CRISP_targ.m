%code to 1) extract endogenous candidate Cas9 binding sites from a target sequence, and 2) screen them for orthogonality to context sequence (ex.: genome) 

%***load sequence context sequence (ex. genome), and target sequence from
%local file
[HeadC, SeqC] = fastaread('E.coli-MG1655-genome.fasta');
%[HeadS, SeqT] = fastaread(File);

%***alternatively load directly form Genbank
%S = getgenbank('M10051');



%generate (reverse) complement
SeqC_C=seqrcomplement(SeqC);
SeqC_RC=seqrcomplement(SeqC);

%convert sequence to vector of integers where A=1, G=3, C=2, T=4
SeqC_v = nt2int(SeqC);
SeqC_C_v=nt2int(SeqC_C);
SeqC_RC_v=nt2int(SeqC_RC);


%***alternative to importing, isolate target sequence from context sequence
%%eventually SeqT_v = getgenbank('gene'); ...

%**define start and length of target sequence
start=4.637613e6; len=716;
%in this case, we are targeting the ArcA gene

%extract target sequence (and its reverse complement) from context sequence
SeqT_v =SeqC_v(start:start+len-1); 
SeqT_C_v =SeqC_C_v(start:start+len-1); 
SeqT_RC_v =fliplr(SeqT_C_v);


%**length of gRNA subsequence 
lseq=23;

%***scan target sequences to find GG positions and output a library of their
%associated gRNA candidates 
[sSeqT_F_v, sSeqT_R_v, sSeqT_F_RC_v, sSeqT_R_RC_v]=GG_library(SeqT_v, SeqT_RC_v, lseq); 

%concatenate the 4 types of gRNA candidates 
sSeqT_glued=[sSeqT_F_v ; sSeqT_R_v; sSeqT_F_RC_v ; sSeqT_R_RC_v ];


%cut out target sequence (minus 22 nt (lseq-1) overhangs) from context sequence
SeqCmT_v_LHS =SeqC_v(1:start+lseq-1); 
SeqCmT_v_RHS =SeqC_v(start+len-(lseq-1):end); 
SeqCmT_RC_v_LHS =fliplr(SeqC_C_v(1:start+lseq-1)); 
SeqCmT_RC_v_RHS =fliplr(SeqC_C_v(start+len-(lseq-1):end)); 

%***with the library of small targetting sequences, scan through context sequence bits, and calculate weighted hamming distance at
%each position 
[m_L, mRC_L, m_R, mRC_R]=scan_SeqC(sSeqT_glued, SeqCmT_v_LHS, SeqCmT_v_RHS, SeqCmT_RC_v_LHS, SeqCmT_RC_v_RHS, lseq); 



