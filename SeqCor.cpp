#include	"FastRegression.cpp"
#include	"dt1sse.cpp"
#include	<iostream>
#include	<fstream>
#include	<unistd.h>
#include	<omp.h>
#include	"t1ha.c"
#include	<map>
vector<RNG>	rng;
vector<float>	yvec;
vector<string>	xstr,	gene;
size_t	dims,	trees,	mtry,	leaf;

bool	load_gene(const	char	*F){
	ifstream	fi(F);
	if(!fi)	return	false;
	string	s;
	for(fi>>s;	!fi.eof();	fi>>s)	gene.push_back(s);
	fi.close();
	return	true;
}

bool	load_data(const	char	*X,	const	char	*Y){
	ifstream	fi(X);
	if(!fi)	return	false;
	string	s;	
	for(fi>>s;	!fi.eof();	fi>>s)	xstr.push_back(s);
	fi.close();
	fi.open(Y);
	if(!fi)	return	false;
	float	y;
	for(fi>>y;	!fi.eof();	fi>>y)	yvec.push_back(y);
	fi.close();
	if(xstr.size()!=yvec.size())	return	false;
	cerr<<X<<'\t'<<Y<<'\t'<<xstr.size()<<'\n';
	return	true;
}

void	gene_correction(void){
	if(gene.size()!=yvec.size()){	cerr<<"signal is NOT adjusted by gene\n";	return;	}
	cerr<<"signal is adjusted by gene\n";
	map<string,	double>	sx,	sn;
	for(size_t	i=0;	i<yvec.size();	i++){	sx[gene[i]]+=yvec[i];	sn[gene[i]]+=1;	}
	for(size_t	i=0;	i<yvec.size();	i++)	yvec[i]-=sx[gene[i]]/sn[gene[i]];
}

void	document(void) {
	cerr<<"Usage:	grna [options] grna input output\n";
	cerr<<"\t-g:	covariant	default=NULL\n";	
	cerr<<"\t-d:	2^dims dimensions	default=9\n";	
	cerr<<"\t-t:	trees	default=896\n";
	cerr<<"\t-m:	mtry	default=2\n";
	cerr<<"\t-l:	leaf	default=cbrt(N)\n";
	exit(0);
}

void	seq2vec(string	&S,	float	*P,	uint64_t	Seed){
	uint64_t	s=t1ha(&Seed,	8,	0);
	for(unsigned	i=0;	i<S.size();	i++){
		uint64_t	h=t1ha(&i,	4,	s);
		for(unsigned	j=1;	j<=3;	j++)	if(i+j<=S.size()){	
			P[t1ha(S.c_str()+i,	j,	s)&(dims-1)]+=1;	
			P[t1ha(S.c_str()+i,	j,	h)&(dims-1)]-=1;
		}
	}
}

int	main(int	ac,	char	**av){
	size_t	t0=time(NULL);
	cerr<<"***********************************\n";
	cerr<<"* SeqCor                          *\n";
	cerr<<"* author: Yi Wang                 *\n";
	cerr<<"* email:  godspeed_china@yeah.net *\n";
	cerr<<"* date:   06/Nov/2018             *\n";
	cerr<<"***********************************\n";
	dims=9;	trees=896;	mtry=2;	leaf=0;
	int	opt;
	while((opt=getopt(ac,	av,	"g:d:t:m:l:"))>=0) {
		switch(opt) {
		case	'g':	load_gene(optarg);	break;
		case	'd':	dims=atoi(optarg);	break;
		case	't':	trees=atoi(optarg);	break;
		case	'm':	mtry=atoi(optarg);	break;
		case	'l':	leaf=atoi(optarg);	break;
		default:	document();
		}
	}
	if(ac<optind+3)	document();
	if(!load_data(av[optind],av[optind+1]))	return	false;
	gene_correction();
	rng.resize(omp_get_num_procs());	rng[0].set(time(NULL));	
	for(size_t	i=1;	i<rng.size();	i++)	rng[i].set(rng[0].get());
	dims=1ULL<<dims;
	cerr<<"dims\t"<<dims<<'\n';
	cerr<<"trees\t"<<trees<<'\n';
	cerr<<"mtry\t"<<mtry<<'\n';
	cerr<<"leaf\t"<<(leaf?leaf:(unsigned)log(xstr.size()))<<'\n';

	vector<float>	ooby(xstr.size());	vector<unsigned>	oobn(xstr.size());	
	#pragma omp parallel for
	for(size_t	t=0;	t<trees;	t++){
		DecisionTree	tree;
		tree.trainn=xstr.size();	tree.testn=0;	tree.feature=dims;	tree.leaf=leaf?leaf:log(tree.trainn);	tree.mtry=mtry;
		tree.resize();
		for(size_t	i=0;	i<xstr.size();	i++)	seq2vec(xstr[i],	tree.xmat+i*dims,	t);
		vector<float>	h(tree.trainn+tree.testn);	vector<bool>	bag;
		tree.estimate(yvec.data(),	h.data(),	bag,	rng[omp_get_thread_num()]);
		for(size_t	i=0;	i<tree.trainn;	i++)	if(!bag[i]){
			#pragma omp atomic
			ooby[i]+=h[i];
			#pragma omp atomic
			oobn[i]++;
		}
		cerr<<'=';
	}
	cerr<<'\n';
	double	sx=0,	sxx=0,	sy=0,	syy=0,	sxy=0,	n=xstr.size();
	FastRegression	fr;	fr.clear();
	for(size_t	i=0;	i<ooby.size();	i++){
		double	h=ooby[i]=oobn[i]?ooby[i]/oobn[i]:0,	y=yvec[i];
		sx+=h;	sxx+=h*h;	sy+=y;	syy+=y*y;	sxy+=h*y;	fr.push2(1,y,h);
	}
	fr.estimate(2);
	ofstream	fo(av[optind+2]);
	for(size_t	i=0;	i<ooby.size();	i++)	fo<<ooby[i]*fr.a+fr.c<<'\n';
	fo.close();
	sx/=n;	sy/=n;
	cout<<av[optind]<<'\t'<<av[optind+1]<<'\t'<<(sxy/n-sx*sy)/sqrt(sxx/n-sx*sx)/sqrt(syy/n-sy*sy)<<'\n';
	cerr<<time(NULL)-t0<<"s\n";
	return	0;
}
