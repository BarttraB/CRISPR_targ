function [sSeqT_F_v,sSeqT_R_v,sSeqT_F_RC_v,sSeqT_R_RC_v]=GG_library(SeqT_v, SeqT_RC_v, lseq); 
%finds  GGs (in  both the top and bottom strand), and their associated (forward and reverse) targeting subsequences 


count=0;
count_RC=0;

%find the GG locations in target sequence
for i1=lseq:length(SeqT_v)-lseq;
    if SeqT_v(i1)==3
        if SeqT_v(i1+1)==3
            count=count+1;
            ind(count)=i1+1;
        end
    end
    if SeqT_RC_v(i1)==3
        if SeqT_RC_v(i1+1)==3
            count_RC=count_RC+1;
            ind_RC(count_RC)=i1+1;
        end
    end
end
    
%preallocate space for forward and reverse targeting subsequences (on top and bottom
%strand)
sSeqT_F_v=zeros(length(ind),lseq);
sSeqT_R_v=zeros(length(ind),lseq);

sSeqT_F_RC_v=zeros(length(ind_RC),lseq);
sSeqT_R_RC_v=zeros(length(ind_RC),lseq);

%make library of these sequences
for i2=1:length(ind);
    sSeqT_F_v(i2,:)=SeqT_v(ind(i2)-(lseq-1):ind(i2));
    sSeqT_R_v(i2,:)=fliplr(SeqT_v(ind(i2)-1:ind(i2)+(lseq-2)));
end

for i3=1:length(ind_RC);
    sSeqT_F_RC_v(i3,:)=SeqT_RC_v(ind_RC(i3)-(lseq-1):ind_RC(i3));
    sSeqT_R_RC_v(i3,:)=fliplr(SeqT_RC_v(ind_RC(i3)-1:ind_RC(i3)+(lseq-2)));
end
