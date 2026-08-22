function [ m_L, m_R ]=scan_SeqC3( sSeqT_glued, lf, SeqCmT_glued_LHS, SeqCmT_glued_RHS, lseq ); 
%scan through context sequence and calculate weighted hamming distance at
%each position (with mostly-excised target sequence) 

%define weight vector (integer values)
for j=1:10
	    jr=j;	
	    w(j)=jr;
end
w(11:lseq-3)=10;
w(lseq-2)=0;
w(lseq-1:lseq)=20;


%print out the longer length library  (useful progress info)
disp(['number of loops for i5 to go through = ' int2str( length(sSeqT_glued(:,lseq)) ) ])


%loop to go through targeting sequence libraries 
for i5=1:length(sSeqT_glued(:,lseq))

    tic %start timer 
    i5 %print out index to show progress
    
    %advance along context sequences loop -LHS
    for i6=1:length(SeqCmT_glued_LHS(1,:))-lseq;
    
        if i<=lf
            sSeq_v_LHS=SeqCmT_glued_LHS(1,i6:i6+lseq-1);
        else
            sSeq_v_LHS=SeqCmT_glued_LHS(2,i6:i6+lseq-1);
        end

        %calculate positions of matches(=1) and mismatches(=0)    
        m=sSeqT_glued(i5,:)==sSeq_v_LHS;

        %scale number of matches by positional weight matrix
        m_L(i5, i6)= dot(w,m);

    end
    
    %advance along context sequences loop -RHS
    for i7=1:length(SeqCmT_glued_RHS(1,:))-lseq;
    
        if i<=lf
            sSeq_v_RHS=SeqCmT_glued_RHS(1,i7:i7+lseq-1);
        else
            sSeq_v_RHS=SeqCmT_glued_RHS(2,i7:i7+lseq-1);
        end
        
        %calculate positions of matches(=1) and mismatches(=0)    
        m=sSeqT_glued(i5,:)==sSeq_v_RHS;

        %scale number of matches by positional weight matrix
        m_R(i5, i7)= dot(w,m);

    end
    
    toc %stop timer 
    %timestamp = datestr(clock, 0) %print date/time
end 


   
