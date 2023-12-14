Y = tsne(dataset);

lda = fitcdiscr(Y,str_labels);
ldaClass = resubPredict(lda);
ldaResubErr = resubLoss(lda);

f = figure;

%gscatter(Y(:,1), Y(:,2), str_labels','rgb','osd');
gscatter(Y(:,1),Y(:,2),str_labels', color, '.', 30);
xlabel('tSNE 1');
ylabel('tSNE 2');

bad = ~strcmp(ldaClass,str_labels');
hold on;
plot(Y(bad,1), Y(bad,2), 'kx');
hold off;

cp = cvpartition(str_labels','KFold',10);
cvlda = crossval(lda,'CVPartition',cp);
ldaCVErr = kfoldLoss(cvlda);

qda = fitcdiscr(Y,str_labels','DiscrimType','quadratic');
qdaResubErr = resubLoss(qda);

cvqda = crossval(qda,'CVPartition',cp);
qdaCVErr = kfoldLoss(cvqda);