#include <windows.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_STEPS 200
#define LINE_SIZE 512
#define STEP_MIN 5.0
#define LOW_THRESHOLD 70.0
#define HIGH_THRESHOLD 180.0

typedef struct {
	double current_glucose;
	double trend;
	double carbs;
	double carb_error;
	double fat;
	double protein;
	double meal_insulin;
	double dia;
	double isf;
	double icr;
	double iob;
	int simulations;
	unsigned seen;
} Input;

enum {
	F_BG=1u<<0, F_TREND=1u<<1, F_CARBS=1u<<2, F_CARBERR=1u<<3,
	F_FAT=1u<<4, F_PROTEIN=1u<<5, F_INSULIN=1u<<6, F_DIA=1u<<7,
	F_ISF=1u<<8, F_ICR=1u<<9, F_IOB=1u<<10
};
#define REQUIRED_MASK ((1u << 11) - 1u)

typedef struct {
	int n;
	double time[MAX_STEPS];
	double central[MAX_STEPS];
	double lower[MAX_STEPS];
	double upper[MAX_STEPS];
	double p_low[MAX_STEPS];
	double p_high[MAX_STEPS];
	double central_min, central_max, envelope_min, envelope_max;
	double t_central_min, t_central_max, t_envelope_min, t_envelope_max;
	double probability_any_low, probability_any_high;
} Result;

unsigned long long rng_state = 1;

void rng_seed( unsigned long s )
{
	rng_state = s ? s : 1;
}

unsigned long long rng_next()
{
	unsigned long long x = rng_state;
	x ^= x >> 12;
	x ^= x << 25;
	x ^= x >> 27;
	rng_state = x;
	return x * 2685821657736338717ULL;
}

double uniform01()
{
	return ((rng_next() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

double normal01()
{
	double u1 = uniform01();
	double u2 = uniform01();
	return sqrt( -2.0 * log( u1 ) ) * cos( 2.0 * M_PI * u2 );
}

double clampd( double x, double lo, double hi )
{
	return x < lo ? lo : (x > hi ? hi : x);
}

char *trim( char* s )
{
	char *e;
	while( isspace( (unsigned char) *s ) )
		s++;
	if ( *s == '\0' )
		return s;
	e = s + strlen( s ) - 1;
	while( e > s && isspace( (unsigned char) *e ) )
		*e-- = '\0';
	return s;
}

int parse_double_strict( const char* s, double* out )
{
	char *end;
	errno = 0;
	*out = strtod( s, &end );
	while( isspace( (unsigned char) *end ) )
		end++;
	return errno == 0 && end != s && *end == '\0' && isfinite( *out );
}

int parse_long_strict( const char* s, long* out )
{
	char* end;
	errno = 0;
	*out = strtol( s, &end, 10 );
	while( isspace( (unsigned char) *end ) )
		end++;
	return errno == 0 && end != s && *end == '\0';
}

void defaults( Input* in )
{
	memset( in, 0, sizeof(*in) );
	in->simulations = 2000;
}

int assign_value( Input* in, const char* key, const char* value, int line )
{
	double d;
	long l;
#define SET_DBL(name, field, flag) if (strcmp(key,name)==0) { if(!parse_double_strict(value,&d)){fprintf(stderr,"行%d: %s は数値で指定してください。\n",line,key);return 0;} in->field=d; in->seen|=flag; return 1; }
	SET_DBL( "current_glucose_mgdl", current_glucose, F_BG )
	SET_DBL( "cgm_trend_mgdl_min", trend, F_TREND )
	SET_DBL( "carbs_g", carbs, F_CARBS )
	SET_DBL( "carb_error_g", carb_error, F_CARBERR )
	SET_DBL( "fat_g", fat, F_FAT )
	SET_DBL( "protein_g", protein, F_PROTEIN )
	SET_DBL( "meal_insulin_units", meal_insulin, F_INSULIN )
	SET_DBL( "insulin_action_minutes", dia, F_DIA )
	SET_DBL( "isf_mgdl_per_unit", isf, F_ISF )
	SET_DBL( "icr_g_per_unit", icr, F_ICR )
	SET_DBL( "iob_units", iob, F_IOB )
#undef SET_DBL
	if ( strcmp( key, "simulations" ) == 0 )
	{
		if ( ! parse_long_strict( value, &l ) )
		{
			fprintf( stderr, "%d行: simulations は整数で指定してください。\n", line );
			return 0;
		}
		in->simulations = (int) l;
		return 1;
	}
	fprintf( stderr, "%d行: 不明なキー '%s'。\n", line, key );
	return 0;
}

int load_input( const char* path, Input* in )
{
	FILE* f = fopen( path, "r" );
	char line[LINE_SIZE];
	int line_no = 0;
	if ( ! f )
	{
		perror( path );
		return 0;
	}
	while( fgets( line, sizeof(line), f ) )
	{
		char *p, *eq, *key, *value;
		line_no++;
		p=trim( line );
		if ( *p == '\0' || *p == '#' )
			continue;
		eq = strchr( p, '=' );
		if ( ! eq )
		{
			fprintf( stderr, "%d行: key=value 形式ではない\n", line_no );
			fclose( f );
			return 0;
		}
		*eq='\0';
		key = trim( p );
		value = trim( eq + 1 );
		if ( ! assign_value( in, key, value, line_no ) )
		{
			fclose( f );
			return 0;
		}
	}
	fclose( f );
	return 1;
}

int in_range( const char* name, double x, double lo, double hi )
{
	if ( x < lo || x > hi )
	{
		fprintf( stderr, "【エラー】%s=%.6g は許容範囲 [%.6g, %.6g] 外\n", name, x, lo, hi );
		return 0;
	}
	return 1;
}

int validate( const Input* in )
{
	int ok = 1;
	if ( (in->seen & REQUIRED_MASK) != REQUIRED_MASK )
	{
		fprintf( stderr, "【エラー】必須パラメーターが設定されていない。\n" );
		return 0;
	}
	// 現在の血糖値 (シミュレーション開始時点の血糖) from リブレ (リブレ表示値)
	ok &= in_range( "current_glucose_mgdl", in->current_glucose, 20, 600 );

	// 血糖の変化率 (現在上昇/下降している速度) from リブレ/CGM履歴 (直近の測定値から計算)
	ok &= in_range( "cgm_trend_mgdl_min", in->trend, -10, 10 );

	// 炭水化物量 (これから摂取する糖質) from 画像解析・栄養表示 (±g)
	ok &= in_range( "carbs_g", in->carbs, 0, 300 );

	// 炭水化物推定誤差 (画像認識結果の不確実性) from 画像解析 (g)
	ok &= in_range( "carb_error_g", in->carb_error, 0, 200 );

	// 脂質 (吸収遅延モデル用) from 画像解析・栄養表示 (g)
	ok &= in_range( "fat_g", in->fat, 0, 300 );

	// タンパク質 (吸収遅延モデル用) from 画像解析・栄養表示 (g)
	ok &= in_range( "protein_g", in->protein, 0, 300 );

	// 食事時インスリン (今回実際に投与すると本人が決めた量) from 本人 (Unit)
	ok &= in_range( "meal_insulin_units", in->meal_insulin, 0, 50 );

	// 作用時間 (インスリン1回の投与が血糖に影響を与えると仮定する作用時間) from 個人の設定 (分)
	ok &= in_range( "insulin_action_minutes", in->dia, 120, 720 );

	// ISF (インスリン感受性係数) from 普段の個人設定 (mg/dL/U)
	ok &= in_range( "isf_mgdl_per_unit", in->isf, 5, 500 );

	// I:C比 (インスリン・炭水化物比) from 普段の個人設定 (g/U)
	ok &= in_range( "icr_g_per_unit", in->icr, 1, 100 );

	// IOB (すでに体内で作用中と推定されるインスリン) from 過去の投与記録 (U) ※前回投与量 × 残り時間 / 作用時間 という線形でも良いが、インスリン作用曲線で計算させるべき
	ok &= in_range( "iob_units", in->iob, 0, 50 );

	// 不確実性を想定したシミュレーション回数 (食事量の推定誤差, 吸収速度のばらつき, インスリン感受性のばらつき, 作用時間のばらつき等々)
	if ( in->simulations < 200 || in->simulations > 20000 )
	{
		fprintf( stderr, "【エラー】シミュレーション範囲 (200～20000) を超えた設定は不可。\n" );
		ok = 0;
	}
	return ok;
}

// ヒューリスティックな累積カーネル (妥当性が検証された生理学的モデルではない)
double weibull_cdf( double t, double tau, double shape )
{
	if ( t <= 0 )
		return 0;
	tau = fmax( tau, 1.0 );
	return 1.0 - exp( -pow( t / tau, shape ) );
}

// Beta(2,4) CDF: 早い段階でピークに達し t=0 でちょうど 0, 期間終了時に 1 となる
double insulin_cdf( double t, double duration )
{
	double x;
	if ( t <= 0 )
		return 0;
	if ( t >= duration )
		return 1;
	x = t / duration;
	return 10 * x * x - 20 * x * x * x + 15 * pow( x, 4 ) - 4 * pow( x, 5 );
}

double iob_cdf( double t, double duration )
{
	double x;
	if ( t <= 0 )
		return 0;
	if ( t >= duration )
		return 1;
	x = t / duration;
	return 1 - (1 - x) * (1 - x);
}

double meal_fraction( double t, double fat, double protein, double absorption_scale )
{
	double delayed = clampd( 0.10 + 0.006 * fat + 0.003 * protein, 0.10, 0.65 );
	double fast_tau = 45.0 * absorption_scale;
	double lag = (35.0 + 0.8 * fat + 0.3 * protein) * absorption_scale;
	double slow_tau = (150.0 + 1.5 * fat + 0.6 * protein) * absorption_scale;
	double fast = weibull_cdf( t, fast_tau, 1.7 );
	double slow = weibull_cdf( t - lag, slow_tau, 1.5 );
	return clampd( (1 - delayed) * fast + delayed * slow, 0, 1 );
}

double trend_effect( double t, double trend, double decay_tau )
{
	return trend * decay_tau * (1 - exp( -t / decay_tau ));
}

double predict_one( const Input* in, double t, double carbs, double isf, double icr, double dia, double absorption_scale, double trend, double trend_tau, double model_slope )
{
	double meal = (carbs / icr) * isf * meal_fraction( t, in->fat, in->protein, absorption_scale );
	double insulin = in->meal_insulin * isf * insulin_cdf( t, dia );
	// IOBには投与量に関する情報が含まれておらず、残存効果モデルは0.75*DIAの間で配分
	insulin += in->iob * isf * iob_cdf( t, 0.75 * dia );
	return in->current_glucose + trend_effect( t, trend, trend_tau ) + meal - insulin + model_slope * t;
}

int cmp_double( const void* a, const void* b )
{
	double x = *(const double*) a;
	double y = *(const double*) b;
	return (x > y) - (x < y);
}

double quantile_sorted( const double* a, int n, double p )
{
	double pos = p * (n - 1);
	double frac;
	int i = (int) floor( pos );
	frac = pos - i;
	return i + 1 < n ? a[i] * (1 - frac) + a[i + 1] * frac : a[i];
}

int simulate( const Input* in, Result* r )
{
	int i, j;
	int nsteps = (int) ceil( fmin( 720.0, in->dia * 1.25 ) / STEP_MIN ) + 1;
	double* vals = NULL;
	unsigned char* any_low = NULL;
	unsigned char* any_high = NULL;
	if ( nsteps > MAX_STEPS )
		return 0;
	vals = (double*) malloc( (size_t) in->simulations * nsteps * sizeof(double) );
	any_low = (unsigned char*) calloc( (size_t) in->simulations, 1 );
	any_high = (unsigned char*) calloc( (size_t) in->simulations, 1 );
	if ( ! vals || ! any_low || ! any_high )
	{
		fprintf( stderr, "メモリ確保に失敗。\n" );
		free( vals );
		free( any_low );
		free( any_high );
		return 0;
	}
	rng_seed( GetTickCount() & 0x7FFFFFFF );
	r->n = nsteps;
	for( i = 0; i < nsteps; i++ )
	{
		double t = i * STEP_MIN;
		r->time[i] = t;
		r->central[i] = predict_one( in, t, in->carbs, in->isf, in->icr, in->dia, 1.0, in->trend, 30.0, 0.0 );
	}
	for( j = 0; j < in->simulations; j++ )
	{
		double carb_delta = in->carb_error * ( uniform01() + uniform01() - 1.0 ); /* triangular */
		double carbs = fmax( 0.0, in->carbs + carb_delta );
		double isf = in->isf * clampd( 1.0 + 0.20 * normal01(), 0.50, 1.60 );
		double icr = in->icr * clampd( 1.0 + 0.12 * normal01(), 0.60, 1.50 );
		double dia = clampd( in->dia * ( 1.0 + 0.15 * normal01() ), 120.0, 720.0 );
		double absorption = exp( 0.28 * normal01() );
		double trend = in->trend + 0.30 * normal01();
		double trend_tau = clampd( 30.0 + 10.0 * normal01(), 12.0, 60.0 );
		double slope = 0.025 * normal01(); /* unmodeled drift, mg/dL/min */
		for( i = 0; i < nsteps; i++ )
		{
			double g = predict_one( in, r->time[i], carbs, isf, icr, dia, absorption, trend, trend_tau, slope );
			vals[(size_t)i * in->simulations + j] = g;
			if ( g < LOW_THRESHOLD )
				any_low[j] = 1;
			if ( g > HIGH_THRESHOLD )
				any_high[j] = 1;
		}
	}
	r->probability_any_low = r->probability_any_high = 0;
	for( j = 0; j < in->simulations; j++ )
	{
		r->probability_any_low += any_low[j];
		r->probability_any_high += any_high[j];
	}
	r->probability_any_low /= in->simulations;
	r->probability_any_high /= in->simulations;
	for( i = 0; i < nsteps; i++ )
	{
		double *row = vals + (size_t)i * in->simulations;
		int lowc = 0, highc = 0;
		for( j = 0; j < in->simulations; j++ )
		{
			if ( row[j] < LOW_THRESHOLD )
				lowc++;
			if ( row[j] > HIGH_THRESHOLD )
				highc++;
		}
		r->p_low[i] = (double)lowc / in->simulations;
		r->p_high[i] = (double)highc / in->simulations;
		qsort( row, (size_t)in->simulations, sizeof(double), cmp_double );
		r->lower[i] = quantile_sorted( row, in->simulations, 0.10 );
		r->upper[i] = quantile_sorted( row, in->simulations, 0.90 );
	}
	r->central_min = r->central_max = r->central[0];
	r->envelope_min = r->lower[0];
	r->envelope_max = r->upper[0];
	r->t_central_min = r->t_central_max = r->t_envelope_min = r->t_envelope_max = 0;
	for( i = 1; i < nsteps; i++ )
	{
		if ( r->central[i] < r->central_min )
		{
			r->central_min = r->central[i];
			r->t_central_min = r->time[i];
		}
		if ( r->central[i] > r->central_max )
		{
			r->central_max = r->central[i];
			r->t_central_max = r->time[i];
		}
		if ( r->lower[i] < r->envelope_min )
		{
			r->envelope_min = r->lower[i];
			r->t_envelope_min = r->time[i];
		}
		if ( r->upper[i] > r->envelope_max )
		{
			r->envelope_max = r->upper[i];
			r->t_envelope_max = r->time[i];
		}
	}
	free( vals );
	free( any_low );
	free( any_high );
	return 1;
}

void write_metadata( FILE* f, const Input* in )
{
	fprintf( f, "# input: glucose=%.2f mg/dL, trend=%.3f mg/dL/min, carbs=%.2f +/- %.2f g, fat=%.2f g, protein=%.2f g\n", in->current_glucose, in->trend, in->carbs, in->carb_error, in->fat,in->protein );
	fprintf( f, "# user-decided insulin=%.3f U, IOB=%.3f U, DIA=%.1f min, ISF=%.2f mg/dL/U, ICR=%.2f g/U\n", in->meal_insulin, in->iob, in->dia, in->isf, in->icr );
	fprintf( f, "# range=p10..p90 Monte Carlo scenario interval; not a calibrated confidence/prediction interval\n" );
	fprintf( f, "# assumptions=heuristic meal absorption + beta insulin action + decaying CGM trend; fat/protein alter timing only\n");
}

int write_csv( const char *path, const Input* in, const Result* r )
{
	FILE* f = fopen( path, "w" );
	int i;
	if ( ! f )
	{
		perror( path );
		return 0;
	}
	write_metadata( f, in );
	fprintf( f, "time_min,central_mgdl,p10_mgdl,p90_mgdl,prob_below_70,prob_above_180\n" );
	for( i = 0; i < r->n; i++ )
		fprintf( f, "%.0f,%.2f,%.2f,%.2f,%.4f,%.4f\n", r->time[i], r->central[i], r->lower[i], r->upper[i], r->p_low[i], r->p_high[i] );
	fclose( f );
	return 1;
}

void report( const Input* in, const Result* r )
{
	printf( "\n--- 1型糖尿病患者向け糖質変化推定シミュレータ ---\n" );
	printf( "【注意】投与量の推奨・決定・最適化は行わない (治療判断に使用してはならない)\n\n" );
	printf( "入力: 現在 %.1f mg/dL, CGM傾向 %+.2f mg/dL/分, 炭水化物 %.1f±%.1f g, 脂質 %.1f g, タンパク質 %.1f g\n", in->current_glucose, in->trend, in->carbs, in->carb_error, in->fat, in->protein );
	printf( "      ユーザー決定済み食事時インスリン %.2f U, 残存インスリン %.2f U, DIA %.0f分, ISF %.1f mg/dL/U, ICR %.1f g/U\n", in->meal_insulin, in->iob, in->dia, in->isf, in->icr );
	printf( "中心推定: 最低 %.1f mg/dL（%.0f分）, 最高 %.1f mg/dL（%.0f分）\n", r->central_min, r->t_central_min, r->central_max, r->t_central_max );
	printf( "10–90%%範囲の極値: 下端 %.1f mg/dL（%.0f分）, 上端 %.1f mg/dL（%.0f分）\n", r->envelope_min, r->t_envelope_min, r->envelope_max, r->t_envelope_max );
	printf( "シナリオ内で一度でも70未満: %.1f%% / 180超: %.1f%%\n", 100 * r->probability_any_low, 100 * r->probability_any_high );
	if ( r->probability_any_low > 0 )
		printf( "警告: 低血糖域に入る可能性を含むシナリオあり\n" );
	if ( r->probability_any_high > 0 )
		printf( "警告: 高血糖域に入る可能性を含むシナリオあり\n" );
	printf( "仮定: 簡略吸収曲線、簡略インスリン作用曲線、減衰するCGM傾向、脂質/タンパク質は吸収時相のみ変更。\n" );
	printf( "未考慮: 運動、ストレス、体調、注射部位、基礎インスリン、吸収変動の全容、CGM遅延など。\n" );
}

int self_test()
{
	int fail = 0;
	double prev = -1;
	int i;
	if ( fabs( insulin_cdf( 0, 300 ) ) > 1e-12 || fabs( insulin_cdf( 300, 300 ) - 1) > 1e-12 )
		fail++;
	for( i = 0; i <= 300; i++ )
	{
		double x = insulin_cdf( i, 300 );
		if ( x + 1e-12 < prev || x < 0 || x > 1 )
			fail++;
		prev = x;
	}
	if ( meal_fraction( 0, 20, 20, 1.0 ) != 0 )
		fail++;
	if ( !( meal_fraction( 300, 40, 30, 1.0 ) > meal_fraction( 60, 40, 30, 1.0 ) ) )
		fail++;
	if ( fail )
	{
		fprintf( stderr, "SELF-TEST FAIL: %d\n", fail );
		return 1;
	}
	printf( "SELF-TEST PASS\n" );
	return 0;
}

void usage(const char *p)
{
	fprintf(stderr,"usage: %s --input FILE [--csv FILE]\n       %s --self-test\n",p,p);
}

int main( int argc,char **argv )
{
	unsigned int uCP = GetConsoleOutputCP();
	SetConsoleOutputCP( 65001 );

	const char* input = NULL;
	const char* csv = "prediction.csv";
	Input in;
	Result r;
	int i;
	if ( argc == 2 && strcmp( argv[1], "--self-test") == 0 )
		return self_test();
	for( i = 1; i < argc; i++ )
	{
		if ( strcmp( argv[i], "--input" ) == 0 && i + 1 < argc )
			input = argv[++i];
		else if ( strcmp(argv[i], "--csv" ) == 0 && i + 1 < argc )
			csv = argv[++i];
		else
		{
			usage( argv[0] );
			return 2;
		}
	}
	if ( ! input )
	{
		usage( argv[0] );
		return 2;
	}
	defaults( &in );
	if ( ! load_input( input, &in ) || ! validate( &in ) )
		return 2;
	if ( ! simulate( &in, &r ) )
	{
		fprintf( stderr, "シミュレーション失敗\n" );
		return 3;
	}
	if ( ! write_csv( csv, &in, &r ) )
		return 4;
	report( &in, &r );
	printf( "\nCSV: %s\n", csv );

	SetConsoleOutputCP( uCP );
	return 0;
}
