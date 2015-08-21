function [m_L, mRC_L, m_R, mRC_R]=scan_SeqC(sSeqT_glued, SeqCmT_v_LHS, SeqCmT_v_RHS, SeqCmT_RC_v_LHS, SeqCmT_RC_v_RHS, lseq); 
%scan through context sequence and calculate weighted hamming distance at
%each position (with mostly-excised target sequence) 

%define weight vector
for j=1:10
	    jr=j;	
	    w(j)=jr/10;
end
w(11:lseq-3)=1;
w(lseq-2)=0;
w(lseq-1:lseq)=2;


%loop to go through small targeting sequence library
for i5=1:length(sSeqT_glued(:,1))

    tic %start timer 
    i5 %print out index to show progress
    
    %advance along context sequences loop -LHS
    for i6=1:length(SeqCmT_v_LHS)-lseq;
    
        sSeq_v_LHS=SeqCmT_v_LHS(i6:i6+lseq-1);
        sSeq_RC_v_LHS=SeqCmT_RC_v_LHS(i6:i6+lseq-1);

        m=sSeqT_glued(i5,:)==sSeq_v_LHS;
        mRC=sSeqT_glued(i5,:)==sSeq_RC_v_LHS;
 
        m_L(i5, i6)= dot(w,m);
        mRC_L(i5, i6)= dot(w,mRC);
    end
    
    %advance along context sequences loop -RHS
    for i7=1:length(SeqCmT_v_RHS)-lseq;
        
        sSeq_v_RHS=SeqCmT_v_RHS(i7:i7+lseq-1);
        sSeq_RC_v_RHS=SeqCmT_v_RHS(i7:i7+lseq-1);
        
        m=sSeqT_glued(i5,:)==sSeq_v_RHS;
        mRC=sSeqT_glued(i5,:)==sSeq_RC_v_RHS;
 
        m_R(i5, i7)= dot(w,m);
        mRC_R(i5, i7)= dot(w,mRC);
    end
    
    toc %stop timer 
end 
