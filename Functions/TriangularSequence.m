function h = TriangularSequence(h)

disp('Creating sequence...');
h.Seq =struct;

% create list of which stim to present on each trial    
h.Seq.condnum = [];
seq=[];
for i = 1:size(h.Settings.oddprob,1)
    ntrial_per_stimtype(i,:) = floor(h.Settings.ntrials(i)*h.Settings.oddprob(i,:));
    thisseq=[];
    for ns = 1:size(ntrial_per_stimtype,2)
        thisseq = [thisseq repmat(ns,1,ntrial_per_stimtype(i,ns))];
    end
    h.Seq.condnum=[h.Seq.condnum i*ones(1,length(thisseq))];
    seq = [seq thisseq];
end
randind = randperm(length(seq));
seq = seq(randind);
h.Seq.condnum = h.Seq.condnum(randind);

for i = 1:size(h.Settings.oddprob,1)
    condind = h.Seq.condnum==i;
    h.Seq.signal(1:3,condind)=h.Settings.oddballvalue{i}(1)*ones(3,sum(condind));
    h.Seq.signal(1,seq==1 & condind) = h.Settings.oddballvalue{i}(2);
    h.Seq.signal(2,seq==2 & condind) = h.Settings.oddballvalue{i}(2);
    h.Seq.signal(3,seq==3 & condind) = h.Settings.oddballvalue{i}(2);
end

h.Seq.blocks=ones(1,length(h.Seq.condnum));
