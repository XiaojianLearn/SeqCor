#include	<algorithm>
#include	<stdint.h>
#include	<cstring>
#include	"rng.cpp"
#include	<cstdlib>
#include	<string>
#include	<cfloat>
#include	<vector>
using	namespace	std;
typedef	vector<unsigned>::iterator	iter;
typedef	float	v4sf	__attribute__	((__vector_size__	(16)));
typedef	int	v4si	__attribute__	((__vector_size__	(16)));
union	XMM {
	v4si	i;
	v4sf	f;
};

class	DecisionTree {
private:
	inline	void	split(float	*Y,	float	*H,	RNG	&R,	iter	TrB,	iter	TrE,	iter	TeB,	iter	TeE,	float	*S,	int	*N);
public:
	DecisionTree() {	xmat=NULL;	}
	~DecisionTree() {	free(xmat);	}
	float	*xmat;
	uint64_t	trainn,	testn,	feature,	leaf,	mtry;
	void	resize(void);
	void	estimate(float	*Y,	float	*H,	vector<bool>	&B,	RNG	&R);
};

void	DecisionTree::resize(void) {
	uint64_t	total=trainn+testn;
	if(posix_memalign((void**)&xmat,	16,	total*feature*4))	return;
	memset(xmat,	0,	total*feature*4);
}


void	DecisionTree::estimate(float	*Y,	float	*H,	vector<bool>	&B,	RNG	&R) {
	vector<unsigned>	trid,	teid;	B.resize(trainn);
	for(uint64_t	i=0;	i<trainn;	i++)	if(B[i]=R(2))	trid.push_back(i);	else	teid.push_back(i);
	for(uint64_t	i=0;	i<testn;	i++)	teid.push_back(trainn+i);	
	float	*s;	int	*n;
	if(posix_memalign((void**)&s,	16,	2*mtry*feature*sizeof(float)))	return;
	if(posix_memalign((void**)&n,	16,	mtry*feature*sizeof(int)))	return;
	split(Y,	H,	R,	trid.begin(),	trid.end(),	teid.begin(),	teid.end(),	s,	n);
	free(s);	free(n);
}

struct	Predictor {
	float	*data;
	uint64_t	feature;
	unsigned	f;
	float	t0;
	bool	operator()(unsigned	I) {	return	t0<=data[I*feature+f];	}
};

void	DecisionTree::split(float	*Y,	float	*H,	RNG	&R,	iter	TrB,	iter	TrE,	iter	TeB,	iter	TeE,	float	*S,	int	*N) {
	uint64_t	n=TrE-TrB;	double	sy=0;
	if(n<=leaf) {
		for(iter	i=TrB;	i!=TrE;	++i)	sy+=Y[*i];
		float	mean=sy/n;	
		for(iter	i=TeB;	i!=TeE;	i++)	H[*i]=mean;
		return;
	}
	memset(S,	0,	2*mtry*feature*sizeof(float));	memset(N,	0,	mtry*feature*sizeof(int));
	float	*T0=S+mtry*feature;
	Predictor	pre;	pre.data=xmat;	pre.feature=feature;	pre.f=pre.t0=0;
	for(uint64_t	m=0;	m<mtry;	m++)	for(uint64_t	i=0;	i<feature;	i++) T0[m*feature+i]=pre.data[*(TrB+R.get()%n)*feature+i];
	XMM	one= {1,1,1,1};
	for(iter	i=TrB;	i!=TrE;	i++) {
		if(i+1!=TrE)	__builtin_prefetch(xmat+*(i+1)*feature);
		float	*x=xmat+*i*feature,	y=Y[*i];
		sy+=y;	v4sf	vy= {y,y,y,y};
		for(size_t	m=0;	m<mtry;	m++) {
			float	*ts=S+m*feature,	*tt0=T0+m*feature;	int	*tn=N+m*feature;
			for(size_t	j=0;	j<feature;	j+=4) {
				XMM	mask;
				mask.f=__builtin_ia32_cmpleps(*(v4sf*)(tt0+j),	*(v4sf*)(x+j));
				*(v4sf*)(ts+j)=__builtin_ia32_addps(__builtin_ia32_andps(vy,	mask.f),	*(v4sf*)(ts+j));
				mask.f=__builtin_ia32_andps(one.f,	mask.f);
				*(v4si*)(tn+j)=__builtin_ia32_paddd128(mask.i,	*(v4si*)(tn+j));
			}
		}
	}
	double	g0=(sy/n)*sy,		max=-FLT_MAX;
	for(uint64_t	m=0;	m<mtry;	m++) {
		float	*ts=S+m*feature,	*tt0=T0+m*feature;	int	*tn=N+m*feature;
		for(uint64_t	j=0;	j<feature;	j++){
			double	g=-g0;
			if(tn[j])	g+=(ts[j]/tn[j])*ts[j];
			if(tn[j]<(int)n)	g+=((sy-ts[j])/(n-tn[j]))*(sy-ts[j]);
			if(g>max) {	max=g;	pre.f=j;	pre.t0=tt0[j];	}
		}
	}
	iter	trit=partition(TrB, TrE, pre),	teit=partition(TeB, TeE, pre);
	if(trit!=TrB&&trit!=TrE) {
		if(teit>TeB)	split(Y,H,R,TrB,trit,TeB,teit,S,N);
		if(teit<TeE)	split(Y,H,R,trit,TrE,teit,TeE,S,N);
	} 
	else {
		float	mean=sy/n;
		for(iter	i=TeB;	i!=TeE;	i++)	H[*i]=mean;
	}
}

