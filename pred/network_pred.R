path_data = args[1]
path_out = args[2]
path_train = 'train.txt'

require(e1071)
func.stand_row = function(row,MAX_log2ptot=6,MAX_log2ntot=6,MAX_log10_tfnei_size=3){
	if(!is.vector(row)) row = unlist(row)
	# standardization, 1e-7 added because some denominators are 0
	row[grep('\\+[^t]',names(row))] = row[grep('\\+[^t]',names(row))]/(row['+total']+1e-7)
	row[grep('\\-[^t]',names(row))] = row[grep('\\-[^t]',names(row))]/(row['-total']+1e-7)
	# standardize '+total' and '-total' columns: log2 + scale to (0,1)
	row['+total'] = log2(row['+total']+1)
	row['+total'] = pmin(row['+total'],MAX_log2ptot) / MAX_log2ptot
	row['-total'] = log2(row['-total']+1)
	row['-total'] = pmin(row['-total'],MAX_log2ntot) / MAX_log2ntot
	row['tfnei_size'] = log10(row['tfnei_size'])
	row['tfnei_size'] = pmin(row['tfnei_size'],MAX_log10_tfnei_size) / MAX_log10_tfnei_size
	return(row)
}

#training
d = read.table2(path_train)
# standardization
d = cbind(d[,1:2], t(apply(d[,-(1:2)],1,func.stand_row)))
kern = 'radial'
nu = 0.5 
gamma = 2**-5
model = svm(d,scale=T,type='one-classification',kernel=kern,nu=nu,gamma=gamma)
# prediction
f = file(path_data,'r')
o = file(path_out,'w')
header = strsplit(readLines(f,n=1),'\t')[[1]]
columns_use = c("+total", "+up10kb", "+up2kb", "+up500b", "+down500b", "+down2kb", "+down10kb", "-total", "-up10kb", "-up2kb", "-up500b", "-down500b", "-down2kb", "-down10kb", "tfnei_size", "coexp2tf", "coexp2tfnei_Q1", "coexp2tfnei_M", "coexp2tfnei_Q3", " coexp2gcatalog_Q1", "coexp2gcatalog_M", "coexp2gcatalog_Q3")
columns_use_indices = match(columns_use, header)
line = readLines(f,n=1)
i = 0
while(length(line)!=0){ 
	if (i%%10000==0)
		print(sprintf('  %d...',i))
	i = i+1
	s = strsplit(line,'\t')[[1]]
	row = as.numeric(s[columns_use_indices])
	names(row) = columns_use
	# standardization
	row_stand = func.stand_row(row)
	# prediction
	res = predict(model, matrix(row_stand,nrow=1))
	if (res){
		writeLines(paste(c(s[1:2],row_stand),collapse='\t'), o)
	}
	# next
	line = readLines(f,n=1)
}
close(f)
close(o)


