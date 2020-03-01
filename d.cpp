#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-result"
#pragma GCC diagnostic ignored "-Wregister"

#include <sstream>
#include <tuple>
#include <vector>
#include <utility>
// #include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <bitset>
// #include <ctime>        // std::time
// #include <time>        // std::time
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <math.h>     
#include <numeric>      // std::accumulate
#include <queue>		//	queue, priority_queue
#include <deque>
//#include <boost/sort/spreadsort/spreadsort.hpp>
#include <exception>
#include <boost/filesystem.hpp>
#include <sys/time.h>
#include <sys/stat.h>
#include <parallel/algorithm>
#include "matplotlibcpp.h"
#include "defines.cpp"

namespace plt = matplotlibcpp;
using namespace std;

string	file_name	=	
							"d_tough_choices.txt";	//	default
							// "f_libraries_of_the_world.txt";	//	default
							// "e_so_many_books.txt";	//	default
// int		should_read_alter_outs	=	0;
// int		should_read_alter_outs	=	1;


string	alter_out_file	=	
							"";
							// "f.out";
							// "e.out";
int		repeat_times			=	1;
float	randomization			=	0.;

								
								random_device rndd; 
								mt19937 mt_rnd(rndd());
								uniform_real_distribution<> zero_to_one_rnd(0.0, 1.0);


								template <class T>
								void SORT_LARGE_FIRST(vector<T> &v) 
									{
									sort( ALL(v), greater<T>());
									}

								template <class T>
								vector<T> Shuffled_Per_Indices( vector<T> &v, vector<int> &ind )
									{
									vector<T> nu(
												min( SIZE(v), SIZE(ind) ) 
												);
									F_all(i, nu)
										{
										nu	=	v[ ind[i] ];
										}
									return nu;
									}

								template <class T, class U>
								vector<T> Split_Pair_0( vector< pair<T,U> > &v )
									{
									vector<T> nu( SIZE(v) );
									F_all(i,v)
										{
										nu[i]	=	v[i].A;
										}
									return	nu;
									}

								template <class T, class U>
								vector<U> Split_Pair_1( vector< pair<T,U> > &v )
									{
									vector<U> nu( SIZE(v) );
									F_all(i,v)
										{
										nu[i]	=	v[i].B;
										}
									return	nu;
									}

								template <class T, class U>
								vector<T> Split_Tuple_0_of_2( vector< tuple<T,U> > &v )
									{
									vector<T> nu( SIZE(v) );
									F_all(i,v)
										{
										nu[i]	=	std::get<0>( v[i] );
										}
									return	nu;
									}

								template <class T, class U>
								vector<U> Split_Tuple_1_of_2( vector< tuple<T,U> > &v )
									{
									vector<U> nu( SIZE(v) );
									F_all(i,v)
										{
										nu[i]	=	std::get<1>( v[i] );
										}
									return	nu;
									}

								template <class T>
								vector<T> Split_VV( vector< vector<T> > &v , int n )
									{
									vector<T> nu( SIZE(v) );
									F_all(i,v)
										{
										nu[i]	=	v[i][n];
										}
									return	nu;
									}

								VI Range( int f, int t )	//	int range [F , T)
									{
									VI rng(t-f);
									F_from(i,f,t)
										{
										rng[i]	=	i;
										}
									return	rng;
									}

								template <class T>
								void SORT_LARGE_FIRST_INDICES_BY_VALUES_FROM_ARRAY( vector<int> &v, vector<T> &by_values )
									{
									vector< pair<T, int> > srt( SIZE(v) );
									F_all(i, v)
										{
										srt[i].A	=	by_values[ v[i] ];
										srt[i].B	=	v[i];
										}
									SORT_LARGE_FIRST( srt );
									v	=	Split_Pair_1(srt);
									}
								template <class T>
								void SORT_SMALL_FIRST_INDICES_BY_VALUES_FROM_ARRAY( vector<int> &v, vector<T> &by_values )
									{
									vector< pair<T, int> > srt( SIZE(v) );
									F_all(i, v)
										{
										srt[i].A	=	by_values[ v[i] ];
										srt[i].B	=	v[i];
										}
									SORT_SMALL_FIRST( srt );
									v	=	Split_Pair_1(srt);
									}


								template <class T>	
								T	V_Median(vector<T> v)
									{
									vector<T> x = v;
									SORT_SMALL_FIRST(x);
									return	x[ SIZE(x)/2 ];
									}


								struct timeval time_start;
								double		TIME_SINCE(struct timeval back_then)
								{
									struct timeval now; gettimeofday(&now, NULL);
									return	((now.tv_sec  - back_then.tv_sec) * 1000000u + 
											now.tv_usec - back_then.tv_usec) / 1.e6;
								}


								template <class T>
								map<T, int> Count_Occurrences_In_V(vector<T> &v) 
									{
									map<int, int> m;
									F_elem(e,v)
										{
										m[ e ]++;
										}
									return m;
									}

								template <class T>
								vector< pair<T, int> > Add_Index(vector<T> &v) 
									{
									vector< pair<T, int> >	nu;
									F_all(i,v)
										{
										nu.PUSH_XX(v[i],i);
										}
									return	nu;
									}

								template <class T, class U>
								vector< pair<T, U> > Combine(vector<T> &a,  vector<U> &b) 
									{
									vector< pair<T, U> >	nu;
									F_all(i,a)
										{
										nu.PUSH_XX(a[i],b[i]);
										}
									return	nu;
									}


								template <class T>
								void out_line_V(vector<T> &v) 
									{
									F_all(i,v)
										{
										std::cout	<<	v[i]	<<	" ";
										}
									std::cout	<<	END;
									}

								template <class T>
								void out_vals_at_indices(vector<T> &v, VI &ind) 
									{
									F_all(i,ind)
										{
										std::cout	<<	v[ ind[i] ]	<<	" ";
										}
									std::cout	<<	END;
									}

								template <class T>
								int count_positive(vector<T> &v) 
									{
									int	x=0;
									F_E_WHERE(v,e, e>0)
										x++;
										}
									return x;
									}

								// template <class T, class U>
								// vector< pair<T, U> > Split_A(vector<T> &a,  vector<U> &b) 
								// 	{
								// 	vector< pair<T, U> >	nu;
								// 	F_all(i,a)
								// 		{
								// 		nu.PUSH_XX(a[i],b[i]);
								// 		}
								// 	return	nu;
								// 	}


								// template <typename T> 
								// T MAX(T a) 
								// 	{
								//     return a;
								// 	}

								// template <typename T, typename ... Args> 
								// T MAX(T a, Args ... args) 
								// 	{
								//     return max(MAX(args...), a);
								// 	}

								// template <typename T> 
								// T MIN(T a) 
								// 	{
								//     return a;
								// 	}

								// template <typename T, typename ... Args> 
								// T MIN(T a, Args ... args) 
								// 	{
								//     return min(MIN(args...), a);
								// 	}


								template <typename T,typename U>                                                   
								std::pair<T,U> operator+(const std::pair<T,U> & l,const std::pair<T,U> & r) {   
									return {l.first+r.first,l.second+r.second};                                    
								}                                                                                  
								template <typename T,typename U>                                                   
								std::pair<T,U> operator-(const std::pair<T,U> & l,const std::pair<T,U> & r) {   
									return {l.first-r.first,l.second-r.second};                                    
								}                                                                                  

								template <class T, class S, class C>
									S& Container(priority_queue<T, S, C>& q) 
										{
										struct HackedQueue : private priority_queue<T, S, C> 
											{
											static S& Container(priority_queue<T, S, C>& q) 
												{
												return q.*&HackedQueue::c;
												}
											};
										return HackedQueue::Container(q);
										}


								ifstream infile;
								ofstream outf;


								bool	 smaller_B__pair_int_double	(const pair<int,double> &a, const pair<int,double> &b)
									{
										return	( a.B < b.B );	//	smaller
									}
								bool	 smaller_B__pair_int_float	(const pair<int,float> &a, const pair<int,float> &b)
									{
										return	( a.B < b.B );	//	smaller
									}
								bool 	smaller_B__pair_int			(const pair<int,int> &a,	const pair<int,int> &b)
									{
										return	( a.B < b.B );	//	smaller
									}


								bool 	larger_B__pair_int_double	(const pair<int,double> &a, const pair<int,double> &b)
									{
										return	( a.B > b.B );	//	larger 
									}
								bool 	larger_B__pair_int_float	(const pair<int,float> &a,	const pair<int,float> &b)
									{
										return	( a.B > b.B );	//	larger 
									}
								bool 	larger_B__pair_int			(const pair<int,int> &a, 	const pair<int,int> &b)
									{
										return	( a.B > b.B );	//	larger 
									}
								bool 	larger_B__pair_P_II_int			(const pair<P_II,int> &a, 	const pair<P_II,int> &b)
									{
										return	( a.B > b.B );	//	larger 
									}


class	upset	//	unordered pseudoset
{
std::size_t C;
std::vector<int> data;
 public:
    upset(std::size_t R, std::size_t C) : C(C), data(R*C) {} // constructor definition 
    int operator()(size_t r, size_t c) const { // member function definition
        return data[r*C+c];
    }

};

								int num_ones_in_binary(int number) {
									int numOnes = 0;
									unsigned int unumber = static_cast<unsigned int>(number);
									// F(i,6) 
									while (unumber != 0)
									{
										numOnes += unumber & 1;
										unumber >>= 1;
									}
									return numOnes;
								}
								int n_ones_in_binary(int n, int number) {
									int numOnes = 0;
									unsigned int unumber = static_cast<unsigned int>(number);
									F(i,n) 
									// while (unumber != 0)
									{
										numOnes += unumber & 1;
										unumber >>= 1;
									}
									return numOnes;
								}

								void	create_folder(string	folder_path)
									{
									struct stat info;
									stat( folder_path.c_str() ,	&info );
									if( S_ISDIR(info.st_mode) )
										{
										// OUT("exists")
										return;
										}

									if	(	boost::filesystem::create_directory( folder_path.c_str() )	)
										{
											OUT( folder_path )
											OUTnl
										}
									else
										{
											OUT_2( "UNABLE TO CREATE A FOLDER!", folder_path )
											// OUT_2( "UNABLE TO CREATE A FOLDER!", folder_path )
											OUTnl
										}
									}


								template<typename T>
								void print_queue(T& q) 
									{
									while(!q.empty()) 
										{
										std::cout << q.top() << " ";
										q.pop();
										}
									std::cout << '\n';
									}

								template<typename T>
								string V_To_STR(vector<T> & v) 
									{
									string	x="";
									F_elem(e,v)
										{
										x	+=	STR(e) + " ";
										}
									return x;
									}

								template<typename T>
								void	out_VP(const T & a)
								{
									OUT( SIZE(a) )
									F_all(i,a)
									{
										OUT( a[i].A, a[i].B );
									}
									return;
								}
								template<typename T>
								void	out_VVP(const T & a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										out_VP(a[i]);
									}
									return;
								}

								template<typename T>
								void	out_V(const T &a)
								{
									OUT( SIZE(a) )
									F_all(i,a)
									{
										OUT( a[i] );
									}
									return;
								}
								template<typename T>
								void	out_VV(const T &a)
								{
									OUT( SIZE(a) )
									F_all(i,a)
									{
										out_V( a[i] );
									}
									return;
								}
								template<typename T>
								void	write_V(const T &a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										FOUT( a[i] );
									}
									return;
								}
								template<typename T>
								void	write_VV(const T &a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										write_V( a[i] );
									}
									return;
								}
								template<typename T>
								void	write_VVV(const T &a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										write_VV( a[i] );
									}
									return;
								}

								VI		read_VI()
								{
									VI	a;
									IN_I(x)
									F(i,x)
									{
										IN_I_PUSH(a);
									}
									return a;
								}


								template<typename T>
								void	write_VP(const T &a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										FOUT_P( a[i] );
									}
									return;
								}
								template<typename T>
								void	write_VVP(const T &a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										write_VP( a[i] );
									}
									return;
								}
								template<typename T>
								void	write_VVVP(const T &a)
								{
									FOUT( SIZE(a) )
									F_all(i,a)
									{
										write_VVP( a[i] );
									}
									return;
								}

								VVI		read_VVI()
								{
									VVI	a;
									IN_I(x)
									F(i,x)
									{
										VI inside = read_VI();
										a.PUSH(inside);
									}
									return a;
								}
								VVVI	read_VVVI()
								{
									VVVI	a;
									IN_I(x)
									F(i,x)
									{
										VVI inside = read_VVI();
										a.PUSH(inside);
									}
									return a;
								}

								VP_II		read_VP_II()
								{
									VP_II	a;
									IN_I(x)
									F(i,x)
									{
										IN_P_PUSH(a);
									}
									return a;
								}
								VVP_II		read_VVP_II()
								{
									VVP_II	a;
									IN_I(x)
									F(i,x)
									{
										VP_II inside = read_VP_II();
										a.PUSH(inside);
									}
									return a;
								}
								VVVP_II	read_VVVP_II()
								{
									VVVP_II	a;
									IN_I(x)
									F(i,x)
									{
										VVP_II inside = read_VVP_II();
										a.PUSH(inside);
									}
									return a;
								}











int	n_books;
int	n_libs;
int	n_days;

VI	b_scores;
VI	b_scores_original;
VVI lib_book;
VVI book_libs;	//[book] -> [libs]
VVI lib_book_initial_sorted;
VI	signup;
VI	rate;
VI	all_books_fit_in_days;
VI	max_books_possible_to_fit;
VI	max_possible_score;
VI	crappy_libs;


VI	should_try_fit_this_book;




int	fit_book(
				VI	&lib_order
			,	VVI	&book_order
			,	int	with_2nd_order	=	0
			)
	{
	// OUT(COL_B+"FITTING..")

	VI		book_sent_first_day	(n_books, 1000000);
	VVP_II	book_sent_by		(n_books);	//	book -> lib , time
	VI		book_sent_well_times(n_books,0);
	VVI		book_sent_well_by	(n_books);	//	book -> lib
	VVI		books_well_sent_by	(n_libs);
	VVI		multi_books_by		(n_libs);	//	lib	->	book_order[lib][ _i__]

	int	day	=	0;

	F_elem(l,lib_order)
		{
		day	+=	signup[l];
		// int	cur_book	=	0;
		
		F_all_elem(bk_i, bk, book_order[l])
			{
			int	tm	=	day
						+ (bk_i / rate[l]);

			if	(
					tm
				<	n_days
				)
				{
				book_sent_well_times[ bk ] ++;
				books_well_sent_by[l].PUSH( bk );
				book_sent_well_by[ bk ].PUSH( l );
				}


			book_sent_first_day[ bk ]	=	min	(
													book_sent_first_day[ bk ]
												,	tm
												);
			book_sent_by[ bk ].PUSH_XX( l , tm );
			}
		}


	VP_II	libs_worst_fit (n_libs, P_II( 0,-1 ));	//	value, iter in book_order[lib][ iter ]
	VP_II	libs_best_unfit(n_libs, P_II( 0,-1 ));

	day	=	0;
	int	most_diff		=	0;
	int	most_diff_lib	=	0;
	int	better_book_i	=	0;
	int	worse_book_i	=	0;


	F_elem(l,lib_order)
		{
		day	+=	signup[l];
		
		F_all_elem(bk_i, bk, book_order[l])
			{
			int	tm	=	day 
						+ (bk_i / rate[l]);


			if	(
					tm
				<	n_days
				)
				{
				if	(
						book_sent_well_times[bk]
					>	1
					)
					{
					multi_books_by[l].PUSH( bk_i );
					}

				if	(
						libs_worst_fit[l].A
					<	b_scores_original[ bk ]
					)
					{
					libs_worst_fit[l].A	=	b_scores_original[ bk ];
					libs_worst_fit[l].B	=	bk_i;
					}
				}
			else
				if	(
							libs_best_unfit[l].A
						>	b_scores_original[ bk ]
					&&
						!	book_sent_well_times[ bk ]
					)
					{
					libs_best_unfit[l].A	=	b_scores_original[ bk ];
					libs_best_unfit[l].B	=	bk_i;
					}

			}

		if	(
				libs_best_unfit[l].A
			&&
				libs_worst_fit[l].A
			)
			IF_MORE_UPDATE_3(
								libs_best_unfit[l].A	-	libs_worst_fit[l].A
							,	libs_best_unfit[l].B
							,	libs_worst_fit[l].B
							,	l
							,	most_diff
							,	better_book_i
							,	worse_book_i
							,	most_diff_lib
							);
		}
	

	int	count_places	=	0;
	F(bk,n_books)
		{
		count_places	+=	max(0,book_sent_well_times[ bk ]	-	1);
		}
	// OUT(COL_Y+"count_places", count_places)




	if	(most_diff)
		{
		OUT	(
				COL_Y	+	"swapping for better book"
			,	book_order[most_diff_lib][ better_book_i ]
			,	book_order[most_diff_lib][ worse_book_i ]
			);
		swap(
				book_order[most_diff_lib][ better_book_i ]
			,	book_order[most_diff_lib][ worse_book_i ]
			);
		int	worse_book	=	book_order[most_diff_lib][ worse_book_i ];
		int	better_book	=	book_order[most_diff_lib][ better_book_i ];
		should_try_fit_this_book[ better_book ]	=	0;
		should_try_fit_this_book[ worse_book  ]	=	b_scores_original[ worse_book  ];
		return 1;
		}


	VI	potential_books = should_try_fit_this_book;
	F_all(bk, book_sent_first_day)
		{
		if	(
					book_sent_first_day[ bk ]
				<	n_days
			||
				book_sent_first_day[ bk ]
				==	1000000
			)
			potential_books[bk]	=	0;
		}


	// VI	count_spare_space;
	int	fixed_one	=	0;
	while(1)
		{
		int	most_promising_book	=	V_MAX_I(potential_books);
		// OUT(book_sent_first_day[ most_promising_book ], most_promising_book,"book_sent_first_day[ most_promising_book ]")
		int	p_sc	=	potential_books[ most_promising_book ];

		// if(most_promising_book == 27794)
		// 	OUT("book", most_promising_book, book_sent_first_day[most_promising_book])

		potential_books[			most_promising_book ] = 0;
		should_try_fit_this_book[	most_promising_book ] = 0;
		if	(!p_sc)
			{
			// OUT("!p_sc")
			return 0;
			}

		int	best_place				=	-1;
		int	best_multibook_i		=	1000000;
		int	best_potential_book_i	=	1000000;
		int	best_author				=	1000000;

		F_elem( e , book_sent_by[ most_promising_book ] )
				{
				int	lb	=	e.A;
				// int	time	=	e.B;

				if	(
					!SIZE(multi_books_by[ lb ])
					)
					continue;


				IF_MORE_UPDATE_2(
									SIZE(multi_books_by[ lb ])
								,	lb
								,	multi_books_by[ lb ][0]
								,	best_place
								,	best_author
								,	best_multibook_i
								);	
				fixed_one	=	1;
				// break;////
				}

		if	(fixed_one)
			{
			best_potential_book_i	=	FIND_I(book_order[best_author] , most_promising_book);

			// OUT	(
			// 		COL_G	+	"swapping"
			// 	,	book_order[best_author][ best_multibook_i ]
			// 	,	book_order[best_author][ best_potential_book_i ]
			// 	);
			swap(
					book_order[best_author][ best_multibook_i ]
				,	book_order[best_author][ best_potential_book_i ]
				);
			break;
			}
		else
		if	(
				!fixed_one
			&&
				with_2nd_order
			)
			{	//	alternatively, try to free some book
			// OUT("xXxXx")

			int	worst_added_book	=	1000000;
			F_elem(l,lib_order)
				{
				F_elem(bk,books_well_sent_by[l])
					{
					MINIMIZE(worst_added_book, b_scores_original[bk])
					}
				}
				
			OUT("p_sc",p_sc,	"worst_added_book score",worst_added_book)
			if	(p_sc	>	worst_added_book)
				{
				VI	eligible_lb( SIZE(book_sent_by[ most_promising_book ]) );
				F_all(i, book_sent_by[ most_promising_book ])
					{
					eligible_lb[ i ]	=	book_sent_by[ most_promising_book ][i].A;
					}

				int	switched	=	0;

				F_elem(elib, eligible_lb)
					{
					// OUT("elib",elib)
					F_elem_where( bk , books_well_sent_by[ elib ] ,
							SIZE(book_sent_by[ bk ])
						>	1
						)
						{	//	try to make bk a multibook
						F_elem_where( bk_sender , book_sent_well_by[ bk ] ,
							SIZE(multi_books_by[ bk_sender ])
							)
							{
							int	m_bk_i	=	multi_books_by[bk_sender][0];
							int	bk_i	=	FIND_I(book_order[ bk_sender ], bk);
							swap(
									book_order[ bk_sender ][  bk_i]
								,	book_order[ bk_sender ][m_bk_i]
								);
							should_try_fit_this_book[ most_promising_book ]	= p_sc;
							OUT(COL_Y_BG+"swapped")
							return 1;
							}
						}
					}
				}
			//	if that failed
			// F_elem(elib, eligible_lb)
			// 	{
			// 	// OUT("elib",elib)
			// 	F_elem_where( bk , books_well_sent_by[ elib ] ,
			// 			SIZE(book_sent_by[ bk ])
			// 		>	1
			// 		)
			// 		{	//	try to make bk a multibook
			// 		F_elem_where( bk_sender , book_sent_well_by[ bk ] ,
			// 			! SIZE(multi_books_by[ bk_sender ])
			// 			)
			// 			{
			// 			int	m_bk_i	=	multi_books_by[bk_sender][0];
			// 			int	bk_i	=	FIND_I(book_order[ bk_sender ], bk);
			// 			swap(
			// 					book_order[ bk_sender ][  bk_i]
			// 				,	book_order[ bk_sender ][m_bk_i]
			// 				);
			// 			should_try_fit_this_book[ most_promising_book ]	= p_sc;
			// 			OUT(COL_Y_BG+"swapped")
			// 			return 1;
			// 			}
			// 		}
			// 	}


			// F_elem(lib_sending_most_promising, eligible_lb)
			// 	{
			// 	F_elem( bk_from_lib_sending_most_promising , books_well_sent_by[ lib_sending_most_promising ] )
			// 		{
			// 		}
			// 	}

			}
		}
	// OUT(COL_G+"fixed_one",  fixed_one)
	return fixed_one;
	}







pair<int, VI> scored	(
				VI	&lib_order
			,	VVI	&book_order
			)
	{
	VI	book_sent_first_day(n_books, 1000000);
	int	day	=	0;
	F_elem(l,lib_order)
		{
		day	+=	signup[l];
		int	cur_book	=	0;
		
		F_elem(bk, book_order[l])
			{
			book_sent_first_day[ bk ]	=	min	(
													book_sent_first_day[ bk ]
												,	day 
													+ (cur_book / rate[l])
												);
			cur_book++;
			}
		}
	
	int	score	=	0;
	int	good	=	0;
	int	bad		=	0;
	int	spared	=	0;
	VI	added_books(n_books, 0);
	F_all(bk, book_sent_first_day)
		{
		if	(
				book_sent_first_day[ bk ]
			<	n_days
			)
			{
			score	+=	b_scores_original[ bk ];
			added_books[ bk ]	=	1;
			good++;
			}
		else
		if	(
				book_sent_first_day[ bk ]
			==	1000000
			)
			spared++;
		else
			{
			bad++;
			// OUT("bad book", bk,		book_sent_first_day[ bk ])
			}
		}

	// OUT( good, COL_B + "good | bad", bad,	spared,	COL_B+"spared" )
	
	return PAIR(score, added_books);
	}


void	pre_fit	(
					VI	 lib_order
				,	VVI	&lib_book
				,	int	optimizing	=	0
				)
	{
	// OUT("optimizing",optimizing)
	VF	book_value;
	CLONE_A_TO_B(b_scores_original, book_value);
	// VI	books_fit_days	=	all_books_fit_in_days;
	VI	books_left (SIZE(lib_order));
	VI	nulled_libs			(n_libs,0);
	VI	nulled_books		(n_books,0);
	VI	nulled_books_in_lib	(n_libs,0);
	VI	books_can_fit( SIZE(lib_order) );


		{
		int	day	=	0;
		F_all_elem(li,l,lib_order)
			{
			day	+=	signup[l];
			books_can_fit[li]	=	max(0, (n_days - day) * rate[l]);
			books_left[li]	=	SIZE(lib_book[l]);
			}
		}

	int		keep_cycling=1;
	float	cycle	=	0;
	while(keep_cycling)
		{
		keep_cycling	=	0;
		F_all_elem(li,l,lib_order)
			{
			if	(
					! nulled_libs[l]
				&&
						books_can_fit[li]
					>=	books_left[li]
				)
				{
				nulled_libs[l]	=	1;
				cycle	+=	0.00001;
				F_elem_where(b,lib_book[l], book_value[b])
					{
					keep_cycling	=	1;
					if	(!nulled_books[b])
						{
						nulled_books[b]	=	1;
						book_value[b]	=	cycle;
						F_elem( nll, book_libs[b] )
							{
							nulled_books_in_lib[ nll ]++;
							}
						}
					}
				}
			}

		int	day;
		day	=	0;
		F_all_elem(l_i,ll,lib_order)
			{
			day	+=	signup[ll];
			if	(
				! nulled_libs[ll]
				)
				{
				books_left[l_i]			=	max(0,  SIZE(lib_book[ll]) - nulled_books_in_lib[ll]  );
				}
			}
		}


	//	TODO: PRESORT lib_book by book_value for betterment of bi
	F_elem(l,lib_order)
		{
		// out_vals_at_indices(book_value, lib_book[l]);
		// if	(optimizing == 10)
		SORT_LARGE_FIRST_INDICES_BY_VALUES_FROM_ARRAY( lib_book[l] , book_value );
		// SHUFFLE(lib_book[l])
		// out_vals_at_indices(book_value, lib_book[l]);
		// OUTnl
		// SORT_LARGE_FIRST_INDICES_BY_VALUES_FROM_ARRAY( lib_book[l] , b_scores_original );
		// SORT_SMALL_FIRST_INDICES_BY_VALUES_FROM_ARRAY( lib_book[l] , book_value );
		}

	// out_line_V(book_value);

	int	day	=	0;
	VVF	all_books, later_minus, later_0;
	F_all_elem(li,l,lib_order)
		{
		day	+=	signup[l];
		
		F_all_elem(bi,be,lib_book[l])
			{
			// VF	x	=	{book_value[be], books_can_fit[li] - bi - nulled_books_in_lib[l], be, l};
			VF	x	=	{
								(float)book_value[be]
							,	(float)( books_can_fit[li] - bi ) 
							// , (float)( 1 ) 
							,	(float)be
							,	(float)l
						};
			if	(x[0] < 1.)
				later_0		.PUSH( x );
			else
				all_books	.PUSH( x );
			}

		lib_book[l].resize(0);
		}

	SORT_LARGE_FIRST(all_books);
	USET<int>	books_included;
	F_elem(bk, all_books)
		{
		if	(
				books_included.count( bk[2] )
			&&	optimizing
			)
			{
			later_minus.PUSH( bk );
			}
		else
			{
			lib_book[ bk[3] ].PUSH( bk[2] );
			books_included.insert( bk[2] );
			}
		}

	//	after all of the good ones add the rest
	F_elem(bke, later_minus)
		{
		lib_book[ bke[3] ].PUSH( bke[2] );
		}

	//	add the sure ones
	SORT_LARGE_FIRST(later_0);
	F_elem(bke2, later_0)
		{
		lib_book[ bke2[3] ].PUSH( bke2[2] );
		}
	}






void	place_fully_included_libs_at_the_latest	(
													VI	&lib_order
												,	int	force	=	0
												)
	{
	if	( ! SIZE(lib_order)
		||	
			(
				file_name=="d_tough_choices.txt"
			&&	!force
			)
		)
		return;

	int	day	= 0;
	F(l, SIZE(lib_order) -1)
		{
		int	lib		= 	lib_order[l];
		day			+=	signup[lib];
		int	full_time	=	all_books_fit_in_days[lib];

		int	nxt_lib	= lib_order[l+1];
		int	next_d	= day;
		F_I_E_FROM(lib_order, l2, lib2 , l+1)
			next_d	+=	signup[lib2];
			if	(
					n_days - next_d 
				>=	full_time
				)
				{
				swap(lib_order[l2],  lib_order[l2 -1] );
				// OUT("place_fully_included_libs_at_the_latest swapped",l2)
				}
			}
		}
	}




map< VI , int > prev_fit_and_calc;

int	fit_and_calc(
					VI	 lib_order
				,	VVI	&lib_book
				,	VI	 added_libs
				,	int	 pre_fit_only	=	0
				,	int	 force			=	0
				)
	{
	APPEND_TO( lib_order , added_libs );
	// OUT(V_To_STR(lib_order));

	place_fully_included_libs_at_the_latest(lib_order);

	VI	in_store	=	{pre_fit_only};
	APPEND_TO(in_store , lib_order)

	// if	(file_name=="d_tough_choices.txt")
	// 	pre_fit_only	=	1;


	if	(
			!force
		&&
			prev_fit_and_calc.count(in_store)
		)
		{
		return	prev_fit_and_calc[in_store];
		}


	lib_book					=	lib_book_initial_sorted;
	should_try_fit_this_book	=	b_scores_original;
	// out_V(lib_order);
	int	x=0;


	pre_fit( lib_order , lib_book , pre_fit_only );

	if	(! pre_fit_only)
		while	(
				fit_book(
							lib_order
						,	lib_book
						// ,	1
						)
				)
				{
				x++;
				}

	// int	b4_2nd_order	=
	// 					scored	(
	// 								lib_order
	// 							,	lib_book
	// 							).A;
	// should_try_fit_this_book	=	b_scores_original;
	// OUT("re")
	// while	(
	// 		fit_book(
	// 					lib_order
	// 				,	lib_book
	// 				,	1
	// 				)
	// 		)
	// 		{
	// 		x++;
	// 		}
	// OUT("-re")



	// should_try_fit_this_book	=	b_scores_original;
	// OUT("re2")
	// while	(
	// 		fit_book(
	// 					lib_order
	// 				,	lib_book
	// 				,	1
	// 				)
	// 		)
	// 		{
	// 		x++;
	// 		}
	// OUT("-re2")
	// OUT(COL_Y+"b4_2nd_order",b4_2nd_order)

	int	score=0;

	score	=	scored	(
							lib_order
						,	lib_book
						).A;
	// OUT(COL_B + "OPTIMIZING", COL_G+STR(MAX_SCORE),	COL_B+STR(TIME_SINCE(time_start))	,	x)

	prev_fit_and_calc[in_store]	=	score;
	return	score;
	}












pair<VI, VD>	deep_calc	(
					int day
				, 	VI	&unused
				,	VI	 lib_order
				,	VVI	&lib_book
				,	int	depth			=	2
				,	int	pre_fit_only	=	0
				)
	{
	VI	real_scores	(n_libs, 0);
	VD	averages	(n_libs, 0.);
	double	pw			=	0.9;
	int		minus_day	=	0;
	int		base_score	=	0;
	// pw	=	0.6;
	// pw			=	0.3;
	// minus_day	=	day;
	pw	=	0.5;
	// pw	=	0.9;
	// double	pw	=	1.2;

	// if	(SIZE(lib_order))
	// 	base_score	=	fit_and_calc(
	// 									lib_order
	// 								,	lib_book
	// 								,	{}
	// 								,	pre_fit_only
	// 								);

	F_all(i,real_scores)
		if	(
				day	
				+	1
				+	signup[i]
			<	n_days
			&&
				unused[i]
			)
		{
		int	di	=		day	
					+	1
					+	signup[i];

		real_scores[i]	=	0;
		int	score	=	fit_and_calc(
										lib_order
									,	lib_book
									,	{i}
									,	pre_fit_only
									);
		MAXIMIZE(
					real_scores[i]
				,	score
				);
		MAXIMIZE(
					averages[i]
				,	
						(float)	(
									score
								-	base_score
								)
					/	pow	(
							(float)
								(
								di
								-	minus_day
								)
							,	pw
							)
				);

		if	(depth<2)
			continue;


		F_all(j,real_scores)
			if	(
					unused[j]
				&&
					i!=j
				&&
						di
						+	signup[j]
					<	n_days
				)
				{
				int	dj	=		di
								+	signup[j]
									;

				int	score	=	fit_and_calc(
												lib_order
											,	lib_book
											,	{i,j}
											,	pre_fit_only
											);
				MAXIMIZE(
							real_scores[i]
						,	score
						);
				MAXIMIZE(
							averages[i]
						,	
							(
							(float)	(
										score
									-	base_score
									)
							/	pow	(
									(float)
										(
										dj
										-	minus_day
										)
									,	pw
									)
							)
						);

				if	(depth==2)
					continue;


				F_all(z,real_scores)
					if	(
							unused[z]
						&&
							i!=z
						&&
							j!=z
						&&
								dj
								+	signup[z]
							<	n_days

						)
						{
						int	dz	=		dj
										+	signup[z]
										;

						int	score	=	fit_and_calc(
														lib_order
													,	lib_book
													,	{i,j,z}
													,	pre_fit_only
													);
						MAXIMIZE(
									real_scores[i]
								,	score
								);
						MAXIMIZE(
									averages[i]
								,	
									(
									(float)	(
												score
											-	base_score
											)
									/	pow	(
											(float)
												(
												dz
												-	minus_day
												)
											,	pw
											)
									)
								);

						if	(depth==3)
							continue;




						F_all(m,real_scores)
							if	(
									unused[m]
								&&
									i!=m
								&&
									j!=m
								&&
									z!=m
								&&
										dz
										+	signup[m]
									<	n_days
								)
								{
								int	dm	=		dz
												+	signup[m]
												;

								int	score	=	fit_and_calc(
																lib_order
															,	lib_book
															,	{i,j,z,m}
															,	pre_fit_only
															);
								MAXIMIZE(
											real_scores[i]
										,	score
										);
								MAXIMIZE(
											averages[i]
										,	
											(
											(float)	(
														score
													-	base_score
													)
											/	pow	(
													(float)
														(
														dm
														-	minus_day
														)
													,	pw
													)
											)
										);
								}
						}
				}
		// OUT(COL_Y_BG+"lib", i,	"produces", real_scores[i], signup[i], "signup",	COL_B+STR(TIME_SINCE(time_start)) )
		}

	return	PAIR(real_scores, averages);
	}






pair<int , VI>	wiggle_and_calc	(
								VI	 lib_order
								)
	{
	int	sz	=	SIZE(lib_order);
	int	f	=	mt_rnd()	%	(sz - 1);
	int	s	=	min(
						sz - 1
					,	f	+	1	+	(int)(mt_rnd() % 20)
					);
	swap( lib_order[f], lib_order[s] );
	return	PAIR(
					fit_and_calc(
									lib_order
								,	lib_book
								,	{}
								,	1
								)
				,	lib_order);
	}




pair<int , VI>	change_random_and_calc	(
										VI	 lib_order
										)
	{
	int	sz	=	SIZE(lib_order);
	int	f	=	mt_rnd()	%	sz;
	USET<int>	current_libs( ALL(lib_order) );

	int	rnd_lib	=	mt_rnd()	%	n_libs;
	while	(
				current_libs.count(rnd_lib) 
			||	crappy_libs[ rnd_lib ]
			)
		{
		rnd_lib	=	mt_rnd()	%	n_libs;
		}

	lib_order[f]	=	rnd_lib;


		{
		int	s	=	min(
							sz - 1
						,	ZERO_PLUS(	f	-	10	+	(int)(mt_rnd() % 20)	)
						);
		swap( lib_order[f], lib_order[s] );
		}


	return	PAIR(
					fit_and_calc(
									lib_order
								,	lib_book
								,	{}
								,	1
								)
				,	lib_order);
	}





void	add_random_lib	(
						VI	 lib_order
						)
	{
	int	sz	=	SIZE(lib_order);
	int	f	=	mt_rnd()	%	sz;
	USET<int>	current_libs( ALL(lib_order) );

	int	rnd_lib	=	mt_rnd()	%	n_libs;
	while	(
				current_libs.count(rnd_lib) 
			||	crappy_libs[ rnd_lib ]
			)
		{
		rnd_lib	=	mt_rnd()	%	n_libs;
		}

	lib_order.PUSH(	rnd_lib );
	}







VI	test_b_books_we_have;
int	test_b_n_different_books	=	0;
USET<int>	test_b_cur_libs;
USET<int>	test_b_one_book_lib;
// USET<int>	test_b_two_book_lib;
USET<int>	test_b_critical_libs;
VI	test_b_is_lib_connected;
VI	test_b_uncritical_books;
// VI	test_b_times_removed_lib;

USET<int>	b_books_related_after_removal(int lib, int send_affected_libs = 1)
	{
	test_b_one_book_lib.erase(lib);
	test_b_critical_libs.erase(lib);

	USET<int>	libz;
	USET<int>	libz_affected;
	USET<int>	books_to_rescan;
	F_E_WHERE( lib_book[lib] , book , test_b_books_we_have[book] == 1 )
		F_E_WHERE( book_libs[book], l , test_b_is_lib_connected[l] )
			libz.insert(l);
			}
		}

	F_E_WHERE( lib_book[lib] , book , test_b_books_we_have[book] == 1 ||  test_b_books_we_have[book] == 2 )
		F_E_WHERE( book_libs[book], l , test_b_is_lib_connected[l] )
			libz_affected.insert(l);
			}
		}

	libz.erase(lib);
	// libz_affected	.erase(lib);

	if	(send_affected_libs)
		F_E( libz_affected , l)
			F_E( lib_book[l] , book )
				books_to_rescan.insert( book );
				}
			}


	F_E( libz , l)
		int	is_few_book_lib	=	test_b_critical_libs.count(l);
		if	( is_few_book_lib )
			continue;

		// int	is_one_book_lib	=	test_b_one_book_lib.count(l);

		int	last_books	=	0;
		F_E_WHERE( lib_book[l] , book , test_b_books_we_have[book] == 1 )
			last_books++;
			}

		if( last_books >= 2 )
			{
			test_b_critical_libs.insert( l );
			F_E( lib_book[l] , book )
				test_b_uncritical_books[ book ] ++;
				// if	(
				// 	// 	test_b_books_we_have[book] >= 2
				// 	// &&	
				// 		send_affected_libs
				// 	)
					// books_to_rescan.insert( book );
				}
			}
	

		if( last_books == 1 )
			{
			test_b_one_book_lib.insert( l );
			}

		// if( last_books == 2 )
		// 	test_b_two_book_lib.insert( l );
		// else
		// 	test_b_two_book_lib.erase( l );
		}

	USET<int>	try_these_libs;
	if	(send_affected_libs)
		{
		F_E( books_to_rescan , bk )
			// OUT("DJASKJLKDA")
			// F_E_WHERE( book_libs[bk] , l ,1)//,  ! test_b_critical_libs.count(l)  &&  test_b_is_lib_connected[l] )
			// F_E_WHERE( book_libs[bk] , l ,  ! test_b_critical_libs.count(l)  && ! test_b_one_book_lib.count(l)  &&  test_b_is_lib_connected[l] )
			// F_E_WHERE( book_libs[bk] , l ,  ! test_b_one_book_lib.count(l)  &&  test_b_is_lib_connected[l] )
			F_E_WHERE( book_libs[bk] , l ,  test_b_is_lib_connected[l] )
			// F_E( book_libs[bk] , l)
				try_these_libs.insert(l);
				// OUT("HERE")
				}
			}
		}

	try_these_libs.erase(lib);

	return	try_these_libs;
	}



void	b_books_related_after_adding(int lib)
	{
	USET<int>	libz;
	F_E_WHERE( lib_book[lib] , book , test_b_books_we_have[book] == 2 )
		F_E_WHERE( book_libs[book], l , test_b_is_lib_connected[l] )
			libz.insert(l);
			}
		}

	libz.erase(lib);


	F_E( libz , l )
		int	is_few_book_lib	=	test_b_critical_libs.count(l);
		int	is_one_book_lib	=	test_b_one_book_lib.count(l);
		if	( is_few_book_lib + is_one_book_lib )
				{
				int	last_books	=	0;
				F_E_WHERE( lib_book[l] , book , test_b_books_we_have[book] == 1 )
					last_books++;
					}

				if	(
						last_books < 2
					&&	is_few_book_lib
					)
					{
					test_b_critical_libs.erase( l );
					F_E( lib_book[l] , book )
						test_b_uncritical_books[ book ] --;
						}
					}

				if	(
						last_books == 0
					&&	is_one_book_lib
					)
					{
					test_b_one_book_lib.erase( l );
					}

				if	(
						last_books == 1
					&&	is_few_book_lib
					&&	! is_one_book_lib
					)
					{
					test_b_one_book_lib.insert( l );
					}

				// if( last_books == 2 )
				// 	test_b_two_book_lib.insert( l );
				// else
				// 	test_b_two_book_lib.erase( l );

				}
		}
	}


void	b_add_lib(int lib)
	{
	if	(test_b_is_lib_connected[lib])
		{
		OUTnl
		OUTnl
		OUT(COL_R_BG + "DAFUQ b_add_lib")
		OUTnl
		OUTnl
		terminate();
		return;
		}

	test_b_cur_libs.insert(lib);
	test_b_is_lib_connected[lib]	=	1;

	F_E( lib_book[lib] , bk )
		if	(
			! test_b_books_we_have[bk]
			)
			test_b_n_different_books++;

		test_b_books_we_have[bk]++;
		}

	b_books_related_after_adding(lib);
	}


pair<int, USET<int> >	b_remove_lib(int lib, int send_affected_libs = 1)
	{
	USET<int> newly_critical_books_to_schedule;

	if	( !test_b_is_lib_connected[lib] )
		{
		OUTnl
		OUTnl
		OUT(COL_R_BG + "DAFUQ b_remove_lib")
		OUTnl
		OUTnl
		terminate();
		return PAIR(-1, newly_critical_books_to_schedule);
		}



	test_b_cur_libs.erase(lib);
	test_b_is_lib_connected[lib]	=	0;

	int	critical	=	0;

	F_E( lib_book[lib] , bk )
		test_b_books_we_have[bk]--;

		if	(
				test_b_books_we_have[bk] == 1
			&&
				! test_b_uncritical_books[bk]
			)
			critical++;

		if	(
			! test_b_books_we_have[bk]
			)
			test_b_n_different_books--;
		}

	if( critical == 0 )
		return PAIR(0, newly_critical_books_to_schedule);

	int	critical_libs_before	=	SIZE(test_b_critical_libs);
	int	one_book_lib_before		=	SIZE(test_b_one_book_lib);
	// int	two_book_lib_before		=	SIZE(test_b_two_book_lib);
	newly_critical_books_to_schedule	=	b_books_related_after_removal(lib, send_affected_libs);

	// critical	=	0;
	// F_E( lib_book[lib] , bk )
	// 	if	(
	// 			test_b_books_we_have[bk] == 1
	// 		// &&
	// 		// 	! test_b_uncritical_books[bk]
	// 		)
	// 		critical++;
	// 	}
	
	// return	SIZE(test_b_critical_libs) - test_b_critical_libs_before;
	// return	critical;
	return	{
				// ZERO_PLUS
				(
					critical	*	20
				+	ZERO_PLUS( SIZE(test_b_one_book_lib)	- one_book_lib_before )		*	40
				-	ZERO_PLUS( SIZE(test_b_critical_libs)	- critical_libs_before )	*	20
				// +	ZERO_PLUS( SIZE(test_b_two_book_lib)	- two_book_lib_before )
				// + ( SIZE(test_b_one_book_lib)	- test_b_test_b_one_book_lib_before )
				// - ( SIZE(test_b_critical_libs)	- test_b_critical_libs_before )
				// 0
				)
			,
				newly_critical_books_to_schedule	
			};
	}


void	REMOVAL_preload_pq_um
				(
					int	cutoff
				,	priority_queue< P_PII_I, vector<P_PII_I>, greater<P_PII_I>>	 &pq
				,	vector<P_II>  &um
				,	int	add_random	=	1
				)
	{
	pq	=	{};
	F(i, n_libs)
		{
		if	(! test_b_is_lib_connected[i])
			continue;

		int	sz = SIZE( lib_book[i] );
		int	a,b;
		a	=	test_b_n_different_books;
		auto	[crit, affected_libs]	=	b_remove_lib(i, 0);
		a	-=	test_b_n_different_books;
		b_add_lib(i);

		int	lazy_rnd	=	0;
		int	lazy_rnd2	=	0;
		if	(add_random)
			{
			lazy_rnd	=	LAZY_RND(2000) ;//+ test_b_times_removed_lib[i];
			lazy_rnd2	=	LAZY_RND(1000) ;//+ test_b_times_removed_lib[i];
			}

		if( a<=cutoff )
			{
			pq.push(PAIR_XX_X	(
									crit	+ lazy_rnd
								,	sz 		+ lazy_rnd2
								,	i
								));
			um[ i ]	=	{
							crit	+ lazy_rnd
						,	sz 		+ lazy_rnd2
						};
			}
		}
	}

void	ADD_preload_pq_um
				(
					int	cutoff
				,	priority_queue< P_PII_I> &pq
				,	int	add_random	=	1
				)
	{
	pq	=	{};
	F_WHERE(n_libs, i, !test_b_is_lib_connected[i])
		int	sz = SIZE( lib_book[i] );
		int	a,b;
		a	=	test_b_n_different_books;
		b_add_lib(i);
		a	=	test_b_n_different_books - a;
		b_remove_lib(i, 0);

		int	lazy_rnd	=	0;
		if	(add_random)
			{
			lazy_rnd	=	LAZY_RND(1000);
			}

		if( a >= cutoff )
			{
			pq.push(PAIR_XX_X	(
									a
								,	sz + lazy_rnd
								,	i
								));
			}
		}
	}


void	REMOVAL_greedy	(
							int	cutoff
						,	priority_queue< P_PII_I, vector<P_PII_I>, greater<P_PII_I>>	 &pq
						,	vector<P_II>  &um
						,	int	max_books
						,	int	until_libs
						)
	{
	OUT( COL_R+"REMOVAL_greedy",cutoff)
	while(
			! pq.empty()		
		&&	
			until_libs < SIZE(test_b_cur_libs)
		)
		{
		auto	pqt	=	pq.top();
		pq.pop();
		int	lib				=	pqt.B;
		int	supposed_crit	=	pqt.A.A;

		if	(
			! test_b_is_lib_connected[lib]
			)
			continue;

		if	(
				pqt.A.B
			>	max_books
			)
			continue;

		if	(
				pqt.A
			!=	um[lib]
			)
			{
			// OUT("-old data, !=	um[lib]")
			continue;
			// break;
			}

		int	book_diff	=	test_b_n_different_books;
		auto	[crit, try_next]	=	b_remove_lib(lib);
		book_diff		-=	test_b_n_different_books;

		if(book_diff>cutoff)
			{
			b_add_lib(lib);
			continue;
			}

		if(crit > supposed_crit)
			{
			b_add_lib(lib);
			pq.push(PAIR_XX_X(
								crit
							,	pqt.A.B
							,	lib
							));

			OUT( COL_B+"crit > supposed_crit", crit, supposed_crit)
			// break;
			continue;
			}

		F_E( try_next , l )	//	UPDATES FOR RELATED LIBS
		// F( l , n_libs ){	//	UPDATES FOR RELATED LIBS
			if( ! test_b_is_lib_connected[l] )
				{
				continue;
				}
			int	in_book_diff		=	test_b_n_different_books;
			auto	[in_crit , x_]	=	b_remove_lib(l, 0);
			in_book_diff			-=	test_b_n_different_books;
			int		in_sz			=	SIZE(lib_book[l]);
			P_II	future_um		=	{in_crit , in_sz};
			if	(
					in_book_diff == 0	
				&&	um[l]	!=	future_um 
				)
				{
				pq.push(PAIR(
								future_um
							,	l
							));
				um[l]	=	future_um;
				// OUT("HELLO", l)
				}
			b_add_lib(l);
			}
		// OUT( "crit", crit, pqt.A.B, "lib books | lib", lib	, COL_B+ STR(TIME_SINCE(time_start))  )

		// test_b_times_removed_lib[lib]++;
		}
	}


void	ADD_greedy	(
						int	cutoff
					,	priority_queue< P_PII_I> &pq
					)
	{
	OUT( COL_G+"ADD_greedy",cutoff)
	while( ! pq.empty() )
		{
		auto	pqt	=	pq.top();
		pq.pop();
		int	lib					=	pqt.B;
		int	supposed_book_diff	=	pqt.A.A;

		if	(test_b_is_lib_connected[lib])
			continue;

		int	book_diff	=	test_b_n_different_books;
		b_add_lib(lib);
		book_diff		=	test_b_n_different_books - book_diff;
		// OUT("lib",lib, book_diff, supposed_book_diff, "was supposed")

		if(book_diff < cutoff)
			{
			b_remove_lib(lib, 0);
			continue;
			}

		if(book_diff < supposed_book_diff)
			{
			b_remove_lib(lib, 0);
			pq.push(PAIR_XX_X(
								book_diff
							,	pqt.A.B
							,	lib
							));
			continue;
			}
		}
	}


void	ADD_threshold	(
							int	cutoff
						,	priority_queue< P_PII_I> &pq_larger_first
						,	int	add_random	=	1
						)
	{
	ADD_preload_pq_um(cutoff,pq_larger_first, add_random);
	ADD_greedy(cutoff,pq_larger_first);
	}

void	REMOVAL_threshold	(
								int	cutoff
							,	priority_queue< P_PII_I, vector<P_PII_I>, greater<P_PII_I>>	 &pq
							,	vector<P_II>  &um
							,	int	max_books	=	1000
							,	int	until_libs	=	0
							)
	{
	REMOVAL_preload_pq_um(	cutoff,pq,um,	until_libs==0 );
	REMOVAL_greedy(			cutoff,pq,um, max_books, until_libs);
	}


void	null_in_preperation()
	{
	test_b_cur_libs			=	{};
	test_b_one_book_lib		=	{};
	test_b_critical_libs	=	{};
	test_b_n_different_books=	0;
	test_b_books_we_have.assign(n_books, 0);
	test_b_uncritical_books.assign(n_books, 0);
	test_b_is_lib_connected.assign(n_libs, 0);
	// test_b_times_removed_lib.assign(n_libs, 0);
	}



int main(int argc, char* argv[])
	{
	srand ( unsigned ( time(0) ) );
	// uniform_int_distribution<> rnd_x(0, x) inclusive
		//	EXAMPLE
			// uniform_int_distribution<> rnd_sz	( 
			// 										min(SIZE(targets)-1, 4) 
			// 									,	max	(
			// 												min(SIZE(targets)-1, 4)  + 1
			// 											,	SIZE(targets)/2
			// 											)
			// 									);
			// target_subset.resize( rnd_sz(mt_rnd) );


	OUTnl
	OUTnl
	if	(argc	==	1)
		{
			//	DEBUG	:	setting	filename here
			//	DEBUG	:	setting	filename here
			//	DEBUG	:	setting	filename here
			// file_name	=	"a_example.txt";
			// file_name	=	"c_incunabula.txt";
			// file_name
			// file_name	=	"b_read_on.txt";

			OUT("Please, specify input file as a parameter")
			// return 0;
		}
	else if
		(argc	>	2)	
		{
		string	second_str	=	(string) argv[2];
		if	(second_str	==	"out")
			{
			file_name		=	argv[1];
			alter_out_file	=	(string) argv[3];
			}
		else
			{
			F_from(i,1,argc)
				{
					string c=
								(string) argv[0]
							+	(string) " "  
							+	(string) argv[i];
					cout<<END<<"++++++"<<END;
					OUT(c)
					int syst=system(c.c_str());
				}
			OUT("All files processed.")
			return 0;
			}
		}
	else
		file_name	=	argv[1];




	OUT_2("file_name="		,file_name)
	OUT_2("alter_out_file="	,alter_out_file)
	OUTnl




	string	folder		=
								(string)"S2"
							+	(string)"__"
							+	(string)"__"
							+	"_static"
							// +	(string)"__rnd" + (string)to_string((int)(rand() % 10000))
							;
	string	folder_path	=	"/ramdisk/" + folder + "/";
	create_folder( folder_path	);
	create_folder( folder_path + "sources/"	);
	string	rnd_str	=	(string)to_string((int)(rand() % 10000));
	string	source_subfolder_rnd	=	"sources/sources"
									+	(string)"__rnd_s-" + rnd_str
									+	"."
									+	file_name
									+	"/";

	create_folder( folder_path + source_subfolder_rnd );
	string	copy_sources_directive;
	// copy_sources_directive	=	"cp -R * " + folder_path + "sources/";
	copy_sources_directive	=	"cp * " + folder_path + source_subfolder_rnd 	+	" 2>/dev/null";
	int syst=system( copy_sources_directive.c_str()	);


	infile	.open	("ins/" + file_name);
	outf	.open	(folder_path + "s-"+ rnd_str +"." + file_name + ".out");

	string	rnd_for_name	=	rnd_str;


	
	TIME(time_start);	// TIME_NOW(time_start_2);
	LL	MAX_SCORE	=	0;
	LL	GLOBAL_MAX	=	0;
	VI	SCORE_chart;
	SET_COUT_FIX_PRECISION;
	SET_COUT_PRECISION(4);
	// SOLUTION START







	IN_VAR(n_books)
	IN_VAR(n_libs)
	IN_VAR(n_days)

	book_libs.resize(n_books);

	F(b, n_books)
		{
		IN_I(sc)
		b_scores.PUSH(sc);
		}

	F(lib, n_libs)
		{
		
		IN_I(n_here)
		IN_I(s_t)
		IN_I(aday)
		signup.PUSH(s_t);
		rate.PUSH(aday);
		all_books_fit_in_days.PUSH( (n_here + aday - 1)/aday );
		max_books_possible_to_fit.PUSH( min(n_here , (n_days - s_t)*aday ) );
		VI	x;
		lib_book.PUSH(x);
		F(bk,n_here)
			{
			IN_I(bk_id)
			lib_book[lib].PUSH(bk_id);

			book_libs[bk_id].PUSH(lib);
			}
		}


	//	SORT BOOKS
	F(i,n_libs)
		{
		VP_II	srt_bk;

		F_all(j,lib_book[i])
			{
			// book_libs[j].PUSH(i);

			int	b	=	lib_book[i][j];
			srt_bk.PUSH_XX	(
								b_scores[b]
							,	b
							);
			}
		
		SORT_LARGE_FIRST(srt_bk);

		F_all(j,lib_book[i])
			{
			int	b	=	srt_bk[j].B;
			lib_book[i][j]	=	b;
			}
		}

	crappy_libs.assign(n_libs, 0);
	max_possible_score.resize(n_libs);
	F(i,n_libs)
		{
		max_possible_score[i]	=	0;
		F(j, max_books_possible_to_fit[i] )
			{
			max_possible_score[i]	+=	b_scores[  lib_book[i][j]  ];
			}
		}

	b_scores_original	=	b_scores;
	lib_book_initial_sorted	=	lib_book;
	int	median_signup	=	V_Median(signup);

	// VI	best_lib_order;
	// F
	VI	lib_order;

	VI	GLOBAL_BEST_lib_order;
	VVI	GLOBAL_BEST_lib_book;

	int	count_xxx	=0;

	infile.close();
	if	(! SIZE(alter_out_file))
			{
			int	min_sz		=	10000000;
			int	max_books	=	0;
			// int	n_shuffle	=	0;
			// F(n_shuffle,100000)

			priority_queue< P_PII_I, vector<P_PII_I>, greater<P_PII_I>>	pq;	//	b_remove_lib().A, n lib elements, lib
			priority_queue< P_PII_I>	pq_larger_first;
			vector<P_II>  um(n_libs, {-1,-1});

			int	last_decrease	=	-1;
			float	avg	=	0.;

			F(n_shuffle,1)
				{
				null_in_preperation();
				// F(i, n_libs)
				// 	{
				// 	b_add_lib(i);
				// 	}
				ADD_threshold(4,pq_larger_first, 0);

				VI	results;

				REMOVAL_threshold(0,pq,um,20,15000);
				// REMOVAL_threshold(1,pq,um,20,15000);
				// REMOVAL_threshold(2,pq,um,20,15000);
				// REMOVAL_threshold(3,pq,um,20,15000);
				// REMOVAL_threshold(4,pq,um,20,15000);

				// F(i,60)
				F(i,500000)
					{
					// // REMOVAL_threshold(0,pq,um);
					// REMOVAL_threshold(1,pq,um);

					// if	(i<4000)
					// 	REMOVAL_threshold(2,pq,um);

					// if	(i<2000)
					// 	REMOVAL_threshold(3,pq,um);

					// if	(i<1000)
					// 	REMOVAL_threshold(4,pq,um);

					// if	(i<500)
					// 	REMOVAL_threshold(5,pq,um);


					// if	(i<1000)
					// 	ADD_threshold(3,pq_larger_first);
					// else
					// 	ADD_threshold(1,pq_larger_first);

					// REMOVAL_threshold(0,pq,um);

					// if(i%10==0)
					// 	{
					// 	ADD_threshold(1,pq_larger_first, 0);
					// 	REMOVAL_threshold(0,pq,um, 100, 15000);
					// 	OUT( "stable_injected")
					// 	}
					int	sz;

					// if	(i<	200)
					// 	{
					// 	// REMOVAL_threshold(1,pq,um);
					// 	// REMOVAL_threshold(2,pq,um);
					// 	// REMOVAL_threshold(3,pq,um);
					// 	// REMOVAL_threshold(4,pq,um);
					// 	// REMOVAL_threshold(5,pq,um);
					// 	// REMOVAL_threshold(6,pq,um);
					// 	// REMOVAL_threshold(7,pq,um);
					// 	// REMOVAL_threshold(8,pq,um);
					// 	REMOVAL_threshold(9,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(9,pq_larger_first);
					// 	// REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	400)
					// 	{
					// 	// REMOVAL_threshold(1,pq,um);
					// 	// REMOVAL_threshold(2,pq,um);
					// 	// REMOVAL_threshold(3,pq,um);
					// 	// REMOVAL_threshold(4,pq,um);
					// 	// REMOVAL_threshold(5,pq,um);
					// 	// REMOVAL_threshold(6,pq,um);
					// 	// REMOVAL_threshold(7,pq,um);
					// 	REMOVAL_threshold(8,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(8,pq_larger_first);
					// 	// REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	600)
					// 	{
					// 	// REMOVAL_threshold(1,pq,um);
					// 	// REMOVAL_threshold(2,pq,um);
					// 	// REMOVAL_threshold(3,pq,um);
					// 	// REMOVAL_threshold(4,pq,um);
					// 	// REMOVAL_threshold(5,pq,um);
					// 	// REMOVAL_threshold(6,pq,um);
					// 	REMOVAL_threshold(7,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(7,pq_larger_first);
					// 	// REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	800)
					// 	{
					// 	// REMOVAL_threshold(1,pq,um);
					// 	// REMOVAL_threshold(2,pq,um);
					// 	// REMOVAL_threshold(3,pq,um);
					// 	// REMOVAL_threshold(4,pq,um);
					// 	// REMOVAL_threshold(5,pq,um);
					// 	REMOVAL_threshold(6,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(6,pq_larger_first);
					// 	// REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	1000)
					// 	{
					// 	// REMOVAL_threshold(1,pq,um);
					// 	// REMOVAL_threshold(2,pq,um);
					// 	// REMOVAL_threshold(3,pq,um);
					// 	// REMOVAL_threshold(4,pq,um);
					// 	REMOVAL_threshold(5,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(5,pq_larger_first);
					// 	// REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	12000000)
					// 	{
					// 	REMOVAL_threshold(1,pq,um);
					// 	REMOVAL_threshold(2,pq,um);
					// 	if	(i%400)
					// 		REMOVAL_threshold(3,pq,um);
					// 	REMOVAL_threshold(4,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(4,pq_larger_first);
					// 	REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	140000000)
					// 	{
					// 	REMOVAL_threshold(1,pq,um);
					// 	REMOVAL_threshold(2,pq,um);
					// 	if	(i%400)
					// 		REMOVAL_threshold(3,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(3,pq_larger_first);
					// 	REMOVAL_threshold(0,pq,um);
					// 	}
					// else
					// if	(i<	200000000)
					// 	{
					// 	REMOVAL_threshold(1,pq,um);
					// 	REMOVAL_threshold(2,pq,um);
					// 	OUT( SIZE(test_b_cur_libs) ,"low point")
					// 	ADD_threshold(2,pq_larger_first);
					// 	REMOVAL_threshold(0,pq,um);
					// 	last_decrease	=	i;
					// 	}
					// else
						{
						REMOVAL_threshold(1,pq,um);
						// if(i%41==0)
						// 	REMOVAL_threshold(2,pq,um);
						ADD_threshold(1,pq_larger_first);
						REMOVAL_threshold(0,pq,um);

						if(i%100==0)
							{
							ADD_threshold(1,pq_larger_first, 0);

							REMOVAL_threshold(0,pq,um, 100, 15000);
							OUT( "stable_injected")
							}

						sz	=	SIZE(test_b_cur_libs);

						if	(
							i-last_decrease	>	10000
							)
							break;

						results.PUSH( sz );

						if	(sz<min_sz	&&	sz<15100)
							lib_order.assign( ALL( test_b_cur_libs ) );

						if	(sz<min_sz)
							{
							min_sz			=	sz;
							last_decrease	=	i;
							}

						MINIMIZE( min_sz , sz )
						}

					avg	+=	(float)SIZE(test_b_cur_libs)	/	20.;
					avg	*=	20./21.;

					OUT( COL_Y + STR( SIZE(test_b_cur_libs) ), "test_b_cur_libs" , min_sz, "min_sz|n_shuffle", n_shuffle, "i " + STR(i),	min_sz , avg, COL_B+ STR(TIME_SINCE(time_start)) )
					}

				// out_line_V(results);



				// OUT( SIZE(test_b_cur_libs), "test_b_cur_libs" , min_sz , n_shuffle)

				int	my_books	=	count_positive( test_b_books_we_have );
				// MAXIMIZE( max_books , my_books )

				OUT("min_sz" ,  min_sz, "max_books" ,  COL_Y +STR(max_books),	my_books,	 n_shuffle	, COL_B+ STR(TIME_SINCE(time_start)))
				// if	(max_books == my_books)
				// 	lib_order.assign( ALL( test_b_cur_libs ) );
				}
			// OUT("min_sz" ,  min_sz, "max_books" ,  COL_Y +STR(max_books))
			OUT("min_sz" ,  min_sz)

			// FOUT(SIZE(lib_order))
			// F_all(i,lib_order)
			// 	{
			// 	int	l	=	lib_order[i];
			// 	// outf<< l << ", ";
			// 	cout<< l << ", ";
			// 	}
			// // outf<< END;
			// cout<< END;


			null_in_preperation();
			F_E(lib_order, lib)
				b_add_lib( lib );
				}
			ADD_threshold(1,pq_larger_first, 0);
			REMOVAL_threshold(0,pq,um, 100, 1300);
			ADD_threshold(1,pq_larger_first, 0);
			REMOVAL_threshold(0,pq,um, 100, 15000);
			int	my_books	=	count_positive( test_b_books_we_have );

			lib_order.assign( ALL( test_b_cur_libs ) );
			OUT( SIZE(lib_order) , "SIZE(lib_order) before ones. Books", my_books)

			REMOVAL_threshold(1,pq,um, 100, 15000);

			my_books	=	count_positive( test_b_books_we_have );
			lib_order.assign( ALL( test_b_cur_libs ) );
			OUT( SIZE(lib_order) , "SIZE(lib_order). Books", my_books )

			REMOVAL_threshold(1,pq,um, 100, 13000);
			VI	sorted_libs( ALL(test_b_cur_libs) );
			
			F_E(lib_order, l)
				if	( ! test_b_cur_libs.count(l) )
					sorted_libs.PUSH(l);
				}

			lib_order	=	sorted_libs;
			}
	else
			{
			infile	.open	("alter_outs/" + alter_out_file);
			
			V_0(lib_order)
			IN_I(n_libs_done)
			int	day	=	0;
			F(l,n_libs_done)
				{
				IN_I(lib_num)
				day	+=	signup[lib_num];
				lib_order.PUSH(lib_num);
				OUT( lib_num
					, "lib | day"
					, day
						, "fit_and_calc"
						,	fit_and_calc(
											lib_order
										,	lib_book
										,	{}
										,	1
										)
						);

				VI	read_books;
				// lib_book[ lib_num ].resize(0);

				set<int> books_read;

				IN_I(scan_n)
				F(i,scan_n)
					{
					IN_I(book)
					read_books.PUSH(book);
					books_read.insert(book);
					}
					
				F_elem(bk , lib_book[ lib_num ])
					{
					if	(! books_read.count(bk))
						read_books.PUSH(bk);
					}
				// lib_book[ lib_num ]	=	read_books;
				}
			}



	place_fully_included_libs_at_the_latest(lib_order);
	
	MAX_SCORE	=	fit_and_calc(
						lib_order
					,	lib_book
					,	{}
					// ,	1 //0
					,	0 //0
					,	1 //force. To produce correct output file.
					);

	OUT(SIZE(lib_order),"lib_orders")
	OUT( COL_Y_BG + COL_B +"MAX_SCORE"	, MAX_SCORE, "-> number of books must be", MAX_SCORE/65, 78600 - MAX_SCORE/65 )



	FOUT(SIZE(lib_order))
	F_all(i,lib_order)
		{
		int	l	=	lib_order[i];
		outf<< l;
		// if	(
		// 		i+1
		// 	<	SIZE(lib_order)
		// 	) 
		outf<< " ";
		outf<< SIZE(lib_book[l]);
		outf<< "\n";
		// F_elem(j, lib_book[l])
		// 	{
		// 	outf<< " "<< j;
		// 	}

		F_all(j, lib_book[l])
			{
			// outf<< " "<< j;
			// outf<< j<< " ";
			outf<< lib_book[l][j];
			if	(
					j+1
				<	SIZE(lib_book[l])
				) 
				outf<< " ";

			}

		outf<< "\n";
		}

	// outf<< "\n";

	// F_all(i,lib_order)
	// 	{
	// 	int	l	=	lib_order[i];
	// 	// outf<<	l;
	// 	F_all(j, lib_book[l])
	// 		{
	// 		// outf<< " "<< j;
	// 		// outf<< j<< " ";
	// 		outf<< lib_book[l][j];
	// 		if	(
	// 				j+1
	// 			<	SIZE(lib_book[l])
	// 			) 
	// 			outf<< " ";

	// 		}
	// 	outf<< "\n";
	// 	}

	// out_V(lib_order);

	// SOLUTION END
	outf.close();
	OUTnl
	OUT( COL_Y_BG + COL_B +"MAX_SCORE"	, MAX_SCORE )
	string	rename_comm	=	"mv " + (folder_path + "s-"+ rnd_str +"." + file_name + ".out") + " " + (folder_path + "s-"+ rnd_str +"." + file_name + ".score-" + STR(MAX_SCORE));
	// OUT(rename_comm)
	system(rename_comm.c_str());
	string	zip_comm	=	"zip -r -q " + folder_path + source_subfolder_rnd + "zipped.zip " + folder_path + source_subfolder_rnd + " && echo 'zip success' || echo 'zip failure'";
	// OUT(zip_comm)
	system(zip_comm.c_str());
	// outf.open	(folder_path + "s-"+ rnd_str +"." + file_name + ".score-" + STR(MAX_SCORE));
	// FOUT(MAX_SCORE)
	// outf.close();


	// RAND_RANGE( rnd_r,2,10 );
	// PARTIAL_RAND_RANGE(rnd_r,2,10,2 );
	// SAVE_BAR_FROM_VI_AS(rnd_r,"chart/asd.png");
	OUT("Storing chart")
	// SCORE_chart.PUSH(350000);
	int	max_chart_size	=	20000;
	// int	max_chart_size	=	200;
	if	(
		SIZE(SCORE_chart)
		>	max_chart_size
		)
		{
		OUT("Downsizing chart from", SIZE(SCORE_chart), "to"	,	max_chart_size)
		SCORE_chart.resize(max_chart_size);
		}
	SAVE_BAR_FROM_VI_AS(
							SCORE_chart
						,	"chart/SCORE_chart_" + file_name + "_" + rnd_for_name +".png"
						);

	
	OUTnl

	OUT( "source_subfolder_rnd",source_subfolder_rnd)
	OUT( COL_Y_BG+COL_B+"TIME"	, COL_B+ STR(TIME_SINCE(time_start)) )
	OUTnl
	OUT("file_name=",file_name)
	OUT( folder_path )
	OUTnl
	OUTnl
	return 0;
	}




//	sudo apt-get install libboost-all-dev  python-matplotlib  python-numpy  python2.7-dev
//	sudo mount -t tmpfs -o size=18192M tmpfs /ramdisk
//	echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc

/*
cd hashcode/

export OMP_NUM_THREADS=4



g++ d.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable d_tough_choices.txt




g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable a_example.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable b_read_on.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable c_incunabula.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable d_tough_choices.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable e_so_many_books.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable f_libraries_of_the_world.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable b_read_on.txt c_incunabula.txt d_tough_choices.txt e_so_many_books.txt f_libraries_of_the_world.txt

g++ s.cpp -Ofast -std=c++17  -I/usr/include/python2.7 -lpython2.7 -lboost_system -lboost_filesystem -o executables/executable &&
./executables/executable e_so_many_books.txt out e.out


*/
