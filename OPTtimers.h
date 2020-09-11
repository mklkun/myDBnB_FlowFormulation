#ifndef OPT_HEADER
#define OPT_HEADER

#include "control_panel.h"

#if( OPT_TIMERS )

class OPTtimers {

public:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	/// constructor of the class
	OPTtimers( void )
	{
		ReSet();
	}

	/// start the timer
	void Start( void )
	{
		if( ! ticking )
		{
#if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 2 ) )
			times( &buff );
			t_u = buff.tms_utime;
			t_s = buff.tms_stime;
#elif( ( OPT_TIMERS == 3 ) || ( OPT_TIMERS == 4 ) )
			t_u = times( &buff );
#elif( OPT_TIMERS == 5 )
			t_u = clock();
#endif

			ticking = 1;
		}
	}

	/// stop the timer
	void Stop( void )
	{
		if( ticking )
		{
			Read( u , s );
			ticking = 0;
		}
	}

	/** Return the elapsed time. If the clock is ticking, return the *total*
	time since the last Start() without stopping the clock; otherwise,
	return the total elapsed time of all the past runs of the clock since
	the last ReSet() [see below]. */

	double Read( void )
	{
		double tu = 0;
		double ts = 0;
		Read( tu , ts );
		return( tu + ts );
	}

	/// As Read( void ) but *adds* user and system time to tu and ts.
	void Read( double &tu , double &ts )
	{
		if( ticking )
		{
#if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 2 ) )
			times( &buff );
			tu += ( double( buff.tms_utime - t_u ) ) / double( CLK_TCK );
			ts += ( double( buff.tms_stime - t_s ) ) / double( CLK_TCK );
#elif( ( OPT_TIMERS == 3 ) || ( OPT_TIMERS == 4 ) )
			tu += ( double( times( &buff ) - t_u ) ) / double( CLK_TCK );
#elif( OPT_TIMERS == 5 )
			tu += ( double( clock() - t_u ) ) / double( CLOCKS_PER_SEC );
#endif
		}
		else { tu += u; ts += s; }
	}

	/// reset the timer
	void ReSet( void )
	{
		u = s = 0; ticking = 0;
	}

private:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	double u;      // elapsed *user* time, in seconds
	double s;      // elapsed *system* time, in seconds
	char ticking;  // 1 if the clock is ticking

#if( ( OPT_TIMERS > 0 ) && ( OPT_TIMERS <= 5 ) )
	clock_t t_u;

#if( OPT_TIMERS <= 4 )
	struct tms buff;

#if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 3 ) )
	clock_t t_s;
#endif
#endif
#endif

};  // end( class OPTtimers );

#endif

#endif