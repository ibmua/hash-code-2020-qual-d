#pragma once

// C++17
// #include <sstream>
// #include <tuple>
// #include <vector>
// #include <utility>
// #include <stdlib.h>     /* srand, rand */
// #include <iostream>
// #include <fstream>
// #include <string>
// #include <map>
// #include <bitset>
// #include <ctime>        // std::time
// #include <algorithm>
// #include <random>
// #include <unordered_map>
// #include <unordered_set>
// #include <set>
// #include <math.h>     
// #include <numeric>      // std::accumulate
// #include <queue>		//	queue, priority_queue
// #include <deque>
// //#include <boost/sort/spreadsort/spreadsort.hpp>
// #include <boost/filesystem.hpp>
// #include <sys/time.h>
// #include <sys/stat.h>
// #include <parallel/algorithm>
// #include "matplotlibcpp.h"

// namespace plt = matplotlibcpp;
// using namespace std;


#define		N_ARGS_SEQ(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,N,...) N
#define		N_ARGS(...) N_ARGS_SEQ(__VA_ARGS__, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
#define		CONCAT(a, b)				a ## b
#define		CONCAT_2(a, b)				CONCAT(a, b)

#define		A               				first
#define		B               				second
#define		ALL(a)							(a).begin(),(a).end()
#define		PAIR							make_pair
#define		PAIR_XX_X(a,b,c)				PAIR(PAIR(a,b),c)
#define		PAIR_X_XX(a,b,c)				PAIR(a,PAIR(b,c))
#define		PAIR_XX_XX(a,b,c,d)				PAIR(PAIR(a,b),PAIR(c,d))

#define		SIZE(x)							(int) x.size()

//#define		BOOST_SORT(x)					boost::sort::spreadsort::spreadsort( ALL(x) );
#define		SORT_SMALL_FIRST(x)				sort( ALL(x) );					//	1 2 3 4
#define		SORT_c(x,c)						sort( ALL(x), c );	//	c(a,b) = a<b  =>  1 2 3 4 ...  [true -> a first, false -> b first]



#define		UMAP							unordered_map
#define		USET							unordered_set
#define		VI								vector<int>
#define		VL								vector<long>
#define		LL								long long
#define		LD								long double
#define		VLL								vector<long long>
#define		VLD								vector<long double>
#define		VF								vector<float>
#define		VD								vector<double>
#define		VS								vector< string >
#define		VB								vector< bool >
#define		VVV(a)							vector< vector< vector< a > > >
#define		VVVV(a)							vector< vector< vector< vector< a > > > >
#define		VVVVV(a)						vector< vector< vector< vector< vector< a > > > > >
#define		VVL								vector< vector<long>	>
#define		VVLD							vector< vector<long double>	>
#define		VVI								vector< vector<int>		>
#define		VVVI							vector< VVI		>
#define		VVVVI							vector< VVVI	>
#define		VVB								vector< VB		>
#define		VVVB							vector< VVB		>
#define		VVVVB							vector< VVVB	>
#define		VVF								vector< vector<float>	>
#define		VVD								vector< vector<double>	>
#define		VVS								vector< vector<string>	>
#define		P_IF							pair<int, float>
#define		P_FI							pair<float,	int>
#define		P_DI							pair<double,	int>
#define		P_ID							pair<int, double>
#define		VP_IF							vector< pair<int,float> >
#define		VP_DI							vector< pair<double,int> >
#define		VP_ID							vector< pair<int,double> >
#define		VP_LDI							vector< pair<LD,int> >
#define		VP_ILD							vector< pair<int,LD> >
#define		VP_DD							vector< pair<double, double> >
#define		P_II							pair<int, int>
#define		P_PII_I							pair<P_II	, int>
#define		P_I_PII							pair<int	, P_II>
#define		P_F_PII							pair<float	, P_II>
#define		P_I_PFI							pair<int	, P_FI>
#define		P_I_PDI							pair<int	, P_DI>
#define		P_PII_PII						pair<P_II, P_II>
#define		P_DD							pair<double, double>
#define		VP_II							vector< pair<int,int> >
#define		VT_III							vector< tuple<int,int,int> >
#define		VP_FI							vector< pair<float,int> >
#define		VP_F_PII						vector< P_F_PII >
#define		VP_I_PII						vector< P_I_PII >
#define		VP_PII_I						vector< P_PII_I >
#define		VVP_FI							vector<  VP_FI	>
#define		VVP_II							vector<  VP_II	>
#define		VVVP_II							vector< VVP_II	>
#define		VVVVP_II						vector<VVVP_II	>
#define		P_LLLL							pair<long long, long long>
#define		P_LL							pair<long, long>
#define		C_SIZE(x)						(sizeof(x) / sizeof(x[0]))
#define		CLEAR(a)						memset(a, 0, sizeof(a));
#define		INF								2000000007
#define		PUSH							push_back
#define		PUSH_XX(x,y)					push_back(PAIR(x,y))
#define		PUSH_X_XX(x,y,z)				push_back(PAIR( x , PAIR(y,z) ))
#define		PUSH_XX_X(x,y,z)				push_back(PAIR( PAIR(x,y) , z ))
#define		UNIQ(v)							unique( ALL(v) );
#define		ONLY_UNIQ(v)					v.erase(	UNIQ(v) , v.end()  );
#define		COPY_APPEND_TO(me,v)			{ me.resize(me.size() + v.size());	copy( ALL(v), me.end() -  v.size()); }
#define		APPEND_TO(me,v)					me.insert(	me.end()	, ALL(v) );
#define		APPEND_INTO_FROM_TO(me,v,f,t)	me.insert(	me.end()	, v.begin() + min(SIZE(v),f), v.begin() + min(SIZE(v),t) );
#define		ADD_BEFORE(me,v)				me.insert(	me.begin()	, ALL(v) );
#define		CLONE_A_TO_B(me,v)				{ v.resize(	me.size() );	copy(	ALL(me), v.begin() ); }
#define		CLONE_VV_A_TO_B(a,b)			{ b.resize(  a.size() );	F_all(i_224dk,b){  b[i_224dk].resize(	a[i_224dk].size() ); copy(	ALL(a[i_224dk]), b[i_224dk].begin() ); } }
#define		REMOVE_TILL(me,i)				{ me.erase( me.begin()		, me.begin() + i ); }
#define		REMOVE_AFTER(me,i)				{ me.erase( me.begin() + i	, me.end()		 ); }
#define		REMOVE_RANGE(me,f,t)			{ me.erase( me.begin() + f	, me.begin() + t ); }
#define		REMOVE_I(from,i)				{ from.erase( from.begin() + i ); }
#define		REMOVE_ITER(me,i)				{ me.erase( i ); }
#define		REMOVE_VALUE(me,v)				{ REMOVE_ITER(me, FIND(me,v) ); }
#define		REMOVE_FROM_SET(st,x)			{ st.erase( st.find(x) , st.end() ); }
#define		REVERSE(myvector)				reverse( ALL(myvector) );
#define		V_MAX_I(v)						distance(v.begin(), max_element( ALL(v) ))
#define		V_MIN_I(v)						distance(v.begin(), min_element( ALL(v) ))
#define		V_0(v)							v.resize(0);

#define		MINIMIZE(a,b)					a=min((a),(b));
#define		MAXIMIZE(a,b)					a=max((a),(b));
#define		UNCONDITIONAL_UPDATE(nu, nu_id, best, best_id)							 														{ best = nu; best_id = nu_id; }
#define		UNCONDITIONAL_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)																{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		UNCONDITIONAL_UPDATE_3(nu, nu_id,	nu_2, nu_3, best, best_id, best_2, best_3)													{ best = nu; best_id = nu_id; best_2 = nu_2; best_3 = nu_3; }
#define		IF_LESS_UPDATE(nu, nu_id, best, best_id)								if(nu < best) 											{ best = nu; best_id = nu_id; }
#define		IF_MORE_UPDATE(nu, nu_id, best, best_id)								if(nu > best) 											{ best = nu; best_id = nu_id; }
#define		IF_LESS_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)				if(nu < best) 											{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		IF_LESS_UPDATE_3(nu, nu_id,	nu_2, nu_3, best, best_id, best_2, best_3)	if(nu < best) 											{ best = nu; best_id = nu_id; best_2 = nu_2; best_3 = nu_3; }
#define		IF_MORE_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)				if(nu > best) 											{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		IF_MORE_UPDATE_3(nu, nu_id,	nu_2, nu_3, best, best_id, best_2, best_3)	if(nu > best) 											{ best = nu; best_id = nu_id; best_2 = nu_2; best_3 = nu_3; }
#define		IF_LESS_EQ_UPDATE(nu, nu_id, best, best_id)								if(	nu < best	||	 (nu	==	best) )					{ best = nu; best_id = nu_id; }
#define		IF_MORE_EQ_UPDATE(nu, nu_id, best, best_id)								if(	nu > best	||	 (nu	==	best) )					{ best = nu; best_id = nu_id; }
#define		IF_LESS_EQ_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)			if(	nu < best	||	 (nu	==	best) )					{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		IF_MORE_EQ_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)			if(	nu > best	||	 (nu	==	best) )					{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		IF_LESS_EQ_RND_UPDATE(nu, nu_id, best, best_id)							if(	nu < best	||	((nu	==	best) && (rand()%2)) )	{ best = nu; best_id = nu_id; }
#define		IF_MORE_EQ_RND_UPDATE(nu, nu_id, best, best_id)							if(	nu > best	||	((nu	==	best) && (rand()%2)) )	{ best = nu; best_id = nu_id; }
#define		IF_LESS_EQ_RND_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)		if(	nu < best	||	((nu	==	best) && (rand()%2)) )	{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		IF_MORE_EQ_RND_UPDATE_2(nu, nu_id,	nu_2, best, best_id, best_2)		if(	nu > best	||	((nu	==	best) && (rand()%2)) )	{ best = nu; best_id = nu_id; best_2 = nu_2; }
#define		IF_THEN_UPDATE  (cond, nu_id,		  best_id)							if(	cond )												{ 			 best_id = nu_id; }
#define		IF_THEN_UPDATE_2(cond, nu_id,	nu_2, best_id, best_2)					if(	cond )												{ 			 best_id = nu_id; best_2 = nu_2; }
#define		IF_THEN_UPDATE_3(cond, nu, nu_id,	nu_2, best, best_id, best_2)		if(	cond )												{ best = nu; best_id = nu_id; best_2 = nu_2; }

#define		F(i,n)							for (int i = 0; i < n; ++i)
#define		F_step(i,t,st)					for (int i = 0; i < t; i+=st )
#define		F_from_step(i,f,t,st)			for (int i = f; i < t; i+=st)
#define		F_all(i,n)						for (int i = 0; i < SIZE(n);	++i)
#define		F_all_elem(i,e,n)				auto e=n[0]; for (int i = 0; i < SIZE(n);	++i, e=n[i])
// #define		F_all_elem(i,e,v)				for (const auto &e : v, int i=)	
#define		F_all_except_last(i,n)			for (int i = 0; i < SIZE(n)-1;	++i)
#define		F_all_rev(i,n)					for (int i = SIZE(n)-1; i >= 0; --i)
#define		F_rev(i,n)						for (int i = n-1; i >= 0; --i)
#define		F_rev_till(i,till,n)			for (int i = n-1; i >= till; --i)
#define		F_from(i,f,t)					for (auto i = f; i < t; ++i)
#define		F_from_till(i,f,t)				for (auto i = f; i != t; f<t? i++ : i--)
#define		F_elem(e,v)						for (const auto &e : v)
#define		F_in(it,v)						for (auto it = v.begin() ; it != v.end(); ++it)

#define		F_I(v,i)						for (int i = 0; i < SIZE(v);	++i)		{
#define		F_I_E(v,i,e)					for (int i = 0; i < SIZE(v);	++i)		{	auto	e=v[i];
#define		F_I_E_FROM(v,i,e,f)				for (int i = f; i < SIZE(v);	++i)		{	auto	e=v[i];
#define		F_I_E_FROM_TO(v,i,e,f,t)		for (int i = f; i < t;			++i)		{	auto	e=v[i];
#define		F_I_E_TO(v,i,e,t)				for (int i = 0; i < t;			++i)		{	auto	e=v[i];
#define		F_E(v,e)						for (const auto &e : v)						{
#define		F_E_WHERE(v,e,w)				for (const auto &e : v)				if(w)	{
#define		F_WHERE(n,i,c)					for (int i = 0; i < n; ++i)									if(c){

#define		F_where(i,n,c)					for (int i = 0; i < n; ++i)									if(c)
#define		F_step_where(i,t,st,c)			for (int i = 0; i < t; i+=st )								if(c)
#define		F_from_step_where(i,f,t,st,c)	for (int i = f; i < t; i+=st)								if(c)
#define		F_all_where(i,n,c)				for (int i = 0; i < SIZE(n);	++i)						if(c)
#define		F_all_elem_where(i,e,n,c)		auto e=n[i];	for (int i = 0; i < SIZE(n);	++i)		if(c)
#define		F_all_except_last_where(i,n,c)	for (int i = 0; i < SIZE(n)-1;	++i)						if(c)
#define		F_all_rev_where(i,n,c)			for (int i = SIZE(n)-1; i >= 0; --i)						if(c)
#define		F_rev_where(i,n,c)				for (int i = n-1; i >= 0; --i)								if(c)
#define		F_rev_till_where(i,till,n,c)	for (int i = n-1; i >= till; --i)							if(c)
#define		F_from_where(i,f,t,c)			for (auto i = f; i < t; ++i)								if(c)
#define		F_from_till_where(i,f,t,c)		for (auto i = f; i != t; f<t? i++ : i--)					if(c)
#define		F_elem_where(e,v,c)				for (const auto &e : v)										if(c)
#define		F_iter_where(it,v,c)			for (auto it = v.begin() ; it != v.end(); ++it)				if(c)

#define		RANGE(rng,f,t)					VI rng(t-f); F_from(i_3483,f,t){ rng[i_3483]	=	i_3483; }

#define		ZERO_PLUS(x)					max(0, x)


#define		LAZY_RND(to)								(	mt_rnd()%to			)	//	0 <= "LAZY" RND < to. May not be completely uniform
#define		RND_ONE_IN(x)								(	mt_rnd()%x == 0		)
#define		RND_CHANCE_IS(x)							(	zero_to_one_rnd(mt_rnd) < x	)
#define		RND_CHANCE									(	zero_to_one_rnd(mt_rnd)		)	//	0. <= RND_CHANCE <= 1.
#define		SHUFFLE(r)									shuffle( ALL(r) , mt_rnd);
#define		F_I_TO_SHUFFLED_RANGE(i,t,r)				RANGE(r,0,t); SHUFFLE(r) F_elem(i,r)  //PARTIAL_RAND_RANGE(r,0,t,1) F_elem(i,r)
// #define		PARTIAL_RAND_RANGE(r,f,t,step)				VI r; random_device rd_322; mt19937 mt_rg(rd_322()); { F_from(i,f,t+1){r.PUSH(i);} F_step(i, SIZE(r), step){shuffle(r.begin()+i , r.begin()+ min(i + step, SIZE(r)) , mt_rg);} }
// #define		F_PARTIAL_RAND_RANGE(i,r,f,t,step)			PARTIAL_RAND_RANGE(r,f,t,step) F_elem(i,r)
// #define		RAND_RANGE(r,f,t)							PARTIAL_RAND_RANGE(r,f,t,t-f)	//	[ f , t ] inclusive  //VI r; { F_from(i,f,t){r.PUSH(i);} random_shuffle(ALL(r)); }
// #define		RND_FLOAT					(	(float)rand() / (float)RAND_MAX				)	//	0. <= x <= 1.
// #define		RND_DOUBLE					(	(double)rand() / (double)RAND_MAX			)
// #define		RND_LD						(	(long double)rand() / (long double)RAND_MAX	)

#define		BIN_SEARCH(me,val)				binary_search( ALL(me), val )
#define		BIN_SEARCH_c(me,val,c)			binary_search( ALL(me), val, c )
#define		UP_BOUND_I(V,val)				upper_bound( ALL(V) ,val )    - V.begin()
#define		LOW_BOUND_I(V,val)				lower_bound( ALL(V) ,val )    - V.begin()
#define		UP_BOUND_I_c(V,val,c)			upper_bound( ALL(V) ,val, c ) - V.begin()
#define		LOW_BOUND_I_c(V,val,c)			lower_bound( ALL(V) ,val, c ) - V.begin()
#define		FIND(me,val)					find( ALL(me), val )
#define		FIND_I(v,val)					FIND(v,val)  -  v.begin() 
#define		THIS_EL_INSIDE(k,m)				m.find( k ) != m.end()
#define		THIS_EL_NOT_INSIDE(k,m)			m.find( k ) == m.end()
#define		IF_THIS_EL_INSIDE(k,m)			if(  m.find( k ) != m.end()  )
#define		IF_THIS_EL_NOT_INSIDE(k,m)		if(  m.find( k ) == m.end()  )

#define		THIS_EL_INSIDE_V(e,v)			find( ALL(v) , e ) != v.end()
#define		THIS_EL_NOT_INSIDE_V(e,v)		find( ALL(v) , e ) == v.end()
#define		IF_THIS_EL_INSIDE_V(e,v)		if(  find( ALL(v) , e ) != v.end()  )
#define		IF_THIS_EL_NOT_INSIDE_V(e,v)	if(  find( ALL(v) , e ) == v.end()  )

#define		BAR_VI_FROM_MAP(vi,m)			VI vi; F(i, (*m.rbegin()).A +1 ){vi.PUSH( m[i] ); }  plt::figure_size(1920*2/3, 1080*2/3);  plt::bar( vi );
#define		BAR_FROM_MAP(m)					VI vi; F(i, (*m.rbegin()).A +1 ){vi.PUSH( m[i] ); }  plt::figure_size(1920*2/3, 1080*2/3);  plt::bar( vi );
#define		SAVE_BAR_FROM_MAP_AS(m,name)	VI vi; F(i, (*m.rbegin()).A +1 ){vi.PUSH( m[i] ); }  plt::figure_size(1920*2/3, 1080*2/3);  plt::bar( vi );	plt::save( name );
#define		BAR_FROM_VI(vi)					plt::figure_size(1920*2/3, 1080*2/3);  plt::bar( vi );
#define		SAVE_BAR_FROM_VI_AS(vi,name)	plt::figure_size(1920*2/3, 1080*2/3);  plt::bar( vi );	plt::save( name );
#define		FILL(me,v)						fill(me.begin(), me.end(), v);
#define		ZERO(me)						FILL(me, 0);
#define		CLOCK_TIME_S(f,t)				(t 		 - f) / (double)CLOCKS_PER_SEC
#define		CLOCK_TIME_SINCE(f)				(clock() - f) / (double)CLOCKS_PER_SEC

#define		TIME(now)						gettimeofday(&now, NULL);
#define		TIME_NOW(now)					struct timeval now; gettimeofday(&now, NULL);


#define		STR(a)			to_string(a)
#define		COL_BLK		((string)"\033[30m" )
#define		COL_R		((string)"\033[31m" )
#define		COL_R_BG	((string)"\033[41m" ) + COL_BLK
#define		COL_G		((string)"\033[32m" )
#define		COL_G_BG	((string)"\033[42m" ) + COL_BLK
#define		COL_B		((string)"\033[34m" )
#define		COL_B_BG	((string)"\033[44m" )
#define		COL_Y		((string)"\033[93m" )
#define		COL_Y_BG	((string)"\033[43m" ) + COL_BLK
#define		COL_END		((string)"\033[0m" )

#define		END								std::endl
#define		OUTnl							{std::cout	<<	COL_END << endl;}	//	new line
#define		OUT_LINE						{std::cout	<<	endl;}				//	new line
#define		OUT_1(x)						{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END	<< endl;}
#define		OUT_2(x,y)						{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<< 	COL_END<< endl;}
#define		OUT_3(x,y,z)					{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END	<<	endl;}
#define		OUT_4(x,y,z,k)					{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END	<<	endl;}
#define		OUT_5(x,y,z,k,l)				{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END	<<	endl;}
#define		OUT_6(x,y,z,k,l,u)				{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END+"	"	<<	u	<<	COL_END	<<	endl;}
#define		OUT_7(x,y,z,k,l,u,p)			{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END+"	"	<<	u	<<	COL_END+"	"	<<	p	<<	COL_END	<<	endl;}
#define		OUT_8(x,y,z,k,l,u,p,h)			{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END+"	"	<<	u	<<	COL_END+"	"	<<	p	<<	COL_END+"	"	<<	h	<<	COL_END	<<	endl;}
#define		OUT_9(x,y,z,k,l,u,p,h,b)		{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END+"	"	<<	u	<<	COL_END+"	"	<<	p	<<	COL_END+"	"	<<	h	<<	COL_END+"	"	<<	b	<<	COL_END	<<	endl;}
#define		OUT_10(x,y,z,k,l,u,p,h,b,m)		{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END+"	"	<<	u	<<	COL_END+"	"	<<	p	<<	COL_END+"	"	<<	h	<<	COL_END+"	"	<<	b	<<	COL_END+"	"	<<	m	<<	COL_END	<<	endl;}
#define		OUT_11(x,y,z,k,l,u,p,h,b,m,w)	{std::cout	<<	COL_END+(string)"	"	<<	x	<<	COL_END+"	"	<<	y	<<	COL_END+"	"	<<	z	<<	COL_END+"	"	<<	k	<<	COL_END+"	"	<<	l	<<	COL_END+"	"	<<	u	<<	COL_END+"	"	<<	p	<<	COL_END+"	"	<<	h	<<	COL_END+"	"	<<	b	<<	COL_END+"	"	<<	m	<<	COL_END+"	"	<<	w	<<	COL_END	<<	endl;}
#define		COUT_1(x)						{std::cout	<<	x	<<	endl;}
#define		COUT_2(x,y)						{std::cout	<<	x	<<	" "	<<	y	<< endl;}
#define		COUT_3(x,y,z)					{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	endl;}
#define		COUT_4(x,y,z,k)					{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	endl;}
#define		COUT_5(x,y,z,k,l)				{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	" "	<<	l	<<	endl;}
#define		COUT_6(x,y,z,k,l,u)				{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	" "	<<	l	<<	" "	<<	u	<<	endl;}
#define		COUT_7(x,y,z,k,l,u,p)			{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	" "	<<	l	<<	" "	<<	u	<<	" "	<<	p	<<	endl;}
#define		COUT_8(x,y,z,k,l,u,p,h)			{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	" "	<<	l	<<	" "	<<	u	<<	" "	<<	p	<<	" "	<<	h	<<	endl;}
#define		COUT_9(x,y,z,k,l,u,p,h,b)		{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	" "	<<	l	<<	" "	<<	u	<<	" "	<<	p	<<	" "	<<	h	<<	" "	<<	b	<<	endl;}
#define		COUT_10(x,y,z,k,l,u,p,h,b,m)	{std::cout	<<	x	<<	" "	<<	y	<<	" "	<<	z	<<	" "	<<	k	<<	" "	<<	l	<<	" "	<<	u	<<	" "	<<	p	<<	" "	<<	h	<<	" "	<<	b	<<	" "	<<	m	<<	endl;}

#define		OUT(...)					CONCAT_2(OUT_ , N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		COUT(...)					CONCAT_2(COUT_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		OUT_map(m)					F_elem(e, m ){ OUT_2( e.A, e.B  )}
#define		OUT_MATRIX(m)				{	F_all(i_22d,m){cout<<i_22d<<"	";	F_all(j_3da, m[i_22d] ){cout<<m[i_22d][j_3da]<<" ";}	OUT_LINE	} }
#define		OUT_P(i)					std::cout<<i.A<<" "<<i.B<<END;

#define		SET_COUT_PRECISION(v)				cout.precision(v);
#define		SET_COUT_FIX_PRECISION				cout.setf( std::ios::fixed, std:: ios::floatfield );
#define		SET_FAKE_INPUT_STORE_TRUE(s,t)		t=cout.rdbuf();	istringstream fake_iss(s);	cin.rdbuf( fake_iss.rdbuf());
#define		SET_FAKE_INPUT(s)					istringstream fake_iss(s);	cin.rdbuf(fake_iss.rdbuf());
// #define		SET_FAKE_INPUT(s)				stringstream fake_iss; fake_iss<<s;	cin.rdbuf(fake_iss.rdbuf());
#define		SET_TRUE_INPUT(input)				cin.rdbuf(input);

#define		CIN_S(i)				string		i;			std::cin>>i;
#define		CIN_P1(i)				P_II		i;			std::cin>>i.A>>i.B;
#define		CIN_P2(i)				P_II		i,j;		std::cin>>i.A>>i.B>>j.A>>j.B;
#define		CIN_LL1(i)				long long	i;			std::cin>>i;
#define		CIN_LL2(i,j)			long long	i,j;		std::cin>>i>>j;
#define		CIN_LL3(i,j,k)			long long	i,j,k;		std::cin>>i>>j>>k;
#define		CIN_LL4(i,j,k,x)		long long	i,j,k,x;	std::cin>>i>>j>>k>>x;
#define		CIN_I1(i)				int			i;			std::cin>>i;
#define		CIN_I2(i,j)				int			i,j;		std::cin>>i>>j;
#define		CIN_I3(i,j,k)			int			i,j,k;		std::cin>>i>>j>>k;
#define		CIN_I4(i,j,k,x)			int			i,j,k,x;	std::cin>>i>>j>>k>>x;
#define		CIN_F1(i)				float		i;			std::cin>>i;
#define		CIN_F2(i,j)				float		i,j;		std::cin>>i>>j;
#define		CIN_F3(i,j,k)			float		i,j,k;		std::cin>>i>>j>>k;
#define		CIN_F4(i,j,k,x)			float		i,j,k,x;	std::cin>>i>>j>>k>>x;
#define		CIN_D1(i)				double		i;			std::cin>>i;
#define		CIN_D2(i,j)				double		i,j;		std::cin>>i>>j;
#define		CIN_D3(i,j,k)			double		i,j,k;		std::cin>>i>>j>>k;
#define		CIN_D4(i,j,k,x)			double		i,j,k,x;	std::cin>>i>>j>>k>>x;


#define		IN_LL_1(i)				long long	i;			infile>>i;
#define		IN_LL_2(i,j)			long long	i,j;		infile>>i>>j;
#define		IN_LL_3(i,j,k)			long long	i,j,k;		infile>>i>>j>>k;
#define		IN_LL_4(i,j,k,x)		long long	i,j,k,x;	infile>>i>>j>>k>>x;
#define		IN_S_1(i)				string		i;			infile>>i;
#define		IN_S_2(i)				string		i,j;		infile>>i>>j;
#define		IN_S_3(i)				string		i,j,k;		infile>>i>>j>>k;
#define		IN_S_4(i)				string		i,j,k,x;	infile>>i>>j>>k>>x;
#define		IN_P_1(i)				P_II		i;			infile>>i.A>>i.B;
#define		IN_P_2(i,j)				P_II		i,j;		infile>>i.A>>i.B>>j.A>>j.B;
#define		IN_I_1(i)				int			i;			infile>>i;
#define		IN_I_2(i,j)				int			i,j;		infile>>i>>j;
#define		IN_I_3(i,j,k)			int			i,j,k;		infile>>i>>j>>k;
#define		IN_I_4(i,j,k,x)			int			i,j,k,x;	infile>>i>>j>>k>>x;
#define		IN_D_1(i)				double		i;			infile>>i;
#define		IN_D_2(i,j)				double		i,j;		infile>>i>>j;
#define		IN_D_3(i,j,k)			double		i,j,k;		infile>>i>>j>>k;
#define		IN_D_4(i,j,k,x)			double		i,j,k,x;	infile>>i>>j>>k>>x;
#define		IN_F_1(i)				float		i;			infile>>i;
#define		IN_F_2(i,j)				float		i,j;		infile>>i>>j;
#define		IN_F_3(i,j,k)			float		i,j,k;		infile>>i>>j>>k;
#define		IN_F_4(i,j,k,x)			float		i,j,k,x;	infile>>i>>j>>k>>x;

#define		IN_1(i)					infile>>i;
#define		IN_2(i,j)				infile>>i>>j;
#define		IN_3(i,j,k)				infile>>i>>j>>k;
#define		IN_4(i,j,k,x)			infile>>i>>j>>k>>x;
#define		IN_5(i,j,k,x,o)			infile>>i>>j>>k>>x>>o;

#define		IN_I_PUSH(a)			{IN_I(g)		a.PUSH(g);}
#define		IN_S_PUSH(a)			{IN_S(g)		a.PUSH(g);}
#define		IN_F_PUSH(a)			{IN_F(g)		a.PUSH(g);}
#define		IN_P_PUSH(a)			{IN_P(g)		a.PUSH(g);}

#define		FOUT_P(i)				outf<<i.A<<" "<<i.B<<END;
#define		FOUT_1(i)				outf<<i<<END;
#define		FOUT_2(i,j)				outf<<i<<" "<<j<<END;
#define		FOUT_3(i,j,k)			outf<<i<<" "<<j<<" "<<k<<END;
#define		FOUT_4(i,j,k,x)			outf<<i<<" "<<j<<" "<<k<<" "<<x<<END;
#define		FOUT_5(i,j,k,x,y)		outf<<i<<" "<<j<<" "<<k<<" "<<x<<" "<<y<<END;
#define		FOUT_6(i,j,k,x,y,z)		outf<<i<<" "<<j<<" "<<k<<" "<<x<<" "<<y<<" "<<z<<END;

#define		FOUT(...)				CONCAT_2(FOUT_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		CIN_I(...)				CONCAT_2(CIN_I, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		CIN_LL(...)				CONCAT_2(CIN_LL, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		CIN_F(...)				CONCAT_2(CIN_F, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		CIN_D(...)				CONCAT_2(CIN_D, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_VAR(...)				CONCAT_2(IN_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_I(...)				CONCAT_2(IN_I_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_LL(...)				CONCAT_2(IN_LL_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_D(...)				CONCAT_2(IN_D_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_F(...)				CONCAT_2(IN_F_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_P(...)				CONCAT_2(IN_P_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		IN_S(...)				CONCAT_2(IN_S_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)

#define		FOUT_map(m)				F_elem(e, m ){ FOUT_2( e.A, e.B  )}

// ARRAYS NEED TO BE SORTED FIRST
#define		UNION_A_B_TO(a,b,to)		set_union		(ALL(a),ALL(b),back_inserter(to));
#define		INTERSECTION_A_B_TO(a,b,to)	set_intersection(ALL(a),ALL(b),back_inserter(to));

#define		SUM_VI(v)				accumulate( v.begin(), v.end(), 0) 
#define		SUM_VF(v)				accumulate( v.begin(), v.end(), 0.)
#define		F_AVG_VI(v)				(float) SUM_VI(v)	/	(float) v.size()
#define		F_AVG_VF(v)				SUM_VF(v)			/	(double) v.size()
#define		MEDIAN(v,m)				auto prtly_srtd = v; nth_element(prtly_srtd.begin(), prtly_srtd.begin() + prtly_srtd.size()/2, prtly_srtd.end()); auto m = v[v.size()/2];

#define		ADD_TO_V(v,by)			{ F_all(i_9dx, v){ v[i_9dx]+=by; } }
#define		SUBTRACT_FROM_V(v,by)	{ F_all(i_9dx, v){ v[i_9dx]-=by; } }
#define		MULTIPLY_V_BY(v,by)		{ F_all(i_9dx, v){ v[i_9dx]*=by; } }
#define		DIVIDE_V_BY(v,by)		{ long double by_rev	= (long double)1./(long double)by;	MULTIPLY_V_BY(v,by_rev) }

#define		ADD_TO_VV(v,by)			{ F_all(j_4dx, v){ ADD_TO_V			(v[j_4dx] , by); } }
#define		SUBTRACT_FROM_VV(v,by)	{ F_all(j_4dx, v){ SUBTRACT_FROM_V	(v[j_4dx] , by); } }
#define		MULTIPLY_VV_BY(v,by)	{ F_all(j_4dx, v){ MULTIPLY_V_BY	(v[j_4dx] , by); } }
#define		DIVIDE_VV_BY(v,by)		{ long double by_rev_x = (long double)1./(long double)by;	MULTIPLY_VV_BY(v,by_rev_x) }

#define		ADD_V(v,b)				{ F(i_94kd, min( SIZE(v) , SIZE(b) ) ){ v[i_94kd]+=b[i_94kd]; } }
#define		SUBTRACT_V(v,b)			{ F(i_94kd, min( SIZE(v) , SIZE(b) ) ){ v[i_94kd]-=b[i_94kd]; } }
#define		MULTIPLY_V(v,b)			{ F(i_94kd, min( SIZE(v) , SIZE(b) ) ){ v[i_94kd]*=b[i_94kd]; } }
#define		DIVIDE_BY_V(v,b)		{ F(i_94kd, min( SIZE(v) , SIZE(b) ) ){ v[i_94kd]/=b[i_94kd]; } }
#define		DIVIDE_BY_V_LD(v,b)		{ F(i_94kd, min( SIZE(v) , SIZE(b) ) ){ v[i_94kd] =(long double)v[i_94kd] / (long double)b[i_94kd]; } }
#define		ADD_VV(v,b)				{ F(i_4857, min( SIZE(v) , SIZE(b) ) ){ ADD_V			(v[i_4857] , b[i_4857]) } }
#define		SUBTRACT_VV(v,b)		{ F(i_4857, min( SIZE(v) , SIZE(b) ) ){ SUBTRACT_V		(v[i_4857] , b[i_4857]) } }
#define		MULTIPLY_VV(v,b)		{ F(i_4857, min( SIZE(v) , SIZE(b) ) ){ MULTIPLY_V		(v[i_4857] , b[i_4857]) } }
#define		DIVIDE_BY_VV(v,b)		{ F(i_4857, min( SIZE(v) , SIZE(b) ) ){ DIVIDE_BY_V		(v[i_4857] , b[i_4857]) } }
#define		DIVIDE_BY_VV_LD(v,b)	{ F(i_4857, min( SIZE(v) , SIZE(b) ) ){ DIVIDE_BY_V_LD	(v[i_4857] , b[i_4857]) } }

#define		NORMALIZE_VF(v)				{ float	 mul = 1./( *max_element(ALL(v)) );	MULTIPLY_ALL_BY(v,mul) }
#define		NORMALIZE_VD(v)				{ double mul = 1./( *max_element(ALL(v)) );	MULTIPLY_ALL_BY(v,mul) }
#define		PERCENTAGES_OF_SUM_VF(v)	{ float	 mul = 1./( SUM_VF(v) );			MULTIPLY_ALL_BY(v,mul) }
#define		PERCENTAGES_OF_SUM_VD(v)	{ double mul = 1./( SUM_VF(v) );			MULTIPLY_ALL_BY(v,mul) }


#define		MAX_2(a,b,c)					max(a		 , b)
#define		MAX_3(a,b,c)					max(a		 , max(b,c))
#define		MAX_4(a,b,c,d)					max(max(a,d) , max(b,c))
#define		MAX_5(a,b,c,d,e)				max(max(a,d) , MAX_3(b,c,e))
#define		MAX_6(a,b,c,d,e,f)				max(max(a,b) , MAX_4(c,d,e,f))
#define		MIN_2(a,b,c)					min(a		 , b)
#define		MIN_3(a,b,c)					min(a		 , min(b,c))
#define		MIN_4(a,b,c,d)					min(min(a,d) , min(b,c))
#define		MIN_5(a,b,c,d,e)				min(min(a,b) , MAX_3(c,d,e))
#define		MIN_6(a,b,c,d,e,f)				min(min(a,b) , MAX_4(c,d,e,f))
#define		MAX(...)						CONCAT_2(MAX_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define		MIN(...)						CONCAT_2(MIN_, N_ARGS(__VA_ARGS__))(__VA_ARGS__)

