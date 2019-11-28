#ifndef	RNG_INCLUDED
#define	RNG_INCLUDED
//xoshiro256**
#include	<immintrin.h>
#include	<stdint.h>
#include	<math.h>
class	RNG{
private:
	union	RNG_union_double{	uint64_t	i;	double	f;	};
	union	RNG_union_single{	uint32_t	i;	float	f;	};
	uint64_t	rotl(const	uint64_t	x,	int	k){	return (x << k) | (x >> (64 - k));  }
	double	uniform_positive(void){	RNG_union_double	u;	do	u.i=(get()&0xfffffffffffffull)|0x3ff0000000000000ull;	while(u.f==1.0);	return	u.f-1.0;	}
public:
	uint64_t	s[4],	curr,	round;
	void	set(uint64_t	x){	for(uint64_t	i=0;	i<4;	i++){	uint64_t	z=(x+=0x9e3779b97f4a7c15);	z=(z^(z>>30))*0xbf58476d1ce4e5b9;	z=(z^(z >> 27))*0x94d049bb133111eb;	s[i]=z^(z >> 31);	}	curr=get();	round=0;	}
	uint64_t	get(void){	uint64_t	r=rotl(s[1]*5,7)*9,	t=s[1]<<17;	s[2]^=s[0];	s[3]^=s[1];	s[1]^=s[2];	s[0]^= s[3];	s[2]^=t;	s[3]=rotl(s[3],45);	return r;	}
	double	uniform(void){	RNG_union_double	u;	u.i=(get()&0xfffffffffffffull)|0x3ff0000000000000ull;	return	u.f-1.0;	}
	float	uniform_single(void){	RNG_union_single	u;	u.i=(get()&0x7ffffful)|0x3f800000ul;	return	u.f-1.0f;	}
	double	normal(void){	double	x,	y,	r2;	do{	x=-1.0+2.0*uniform_positive();	y=-1.0+2.0*uniform_positive();	r2=x*x+y*y;	}while(r2>1.0||r2==0);	return	y*sqrt(-2.0*log(r2)/r2);	}
	float	fast_normal(void){	float	s=-4.5f;	RNG_union_single	u;	u.i=(get()&0x7ffffful)|0x3f800000ul;	s+=u.f;	u.i=(get()&0x7ffffful)|0x3f800000ul;	s+=u.f;	u.i=(get()&0x7ffffful)|0x3f800000ul;	s+=u.f;	return	2.0f*s;	}
	bool	operator()(int	K){	
		if(K>0){
			if(round+K>64){	curr=get();	round=0;	}	uint64_t	r=curr&((1ULL<<K)-1ULL);	curr>>=K;	round+=K;	return	r;	
		}
		else{
			K=-K;
			if(round+K>64){	curr=get();	round=0;	}	uint64_t	r=curr&((1ULL<<K)-1ULL);	curr>>=K;	round+=K;	return	!r;	
		}
	}
};
#endif

