/*  Simple PPM functions.  See sample main function for usage.
 *      
 *  int  get_ppm( image * x ) -- reads and allocates PPM ("P6") from stdin.
 *  int  copy_ppm( image x, image * y ) -- allocates a copy of x, puts in *y.
 *
 *  Allocation functions return a nonzero value if successful.
 *
 *  void dealloc_ppm( image x ) -- frees an allocated image.
 *  void put_ppm( image x ) -- writes x in PPM format on stdout.
 *
 *  pixel get_pixel( image z, int x, int y ) --- get RGB pixel at (x,y)
 *  void  put_pixel( image z, int x, int y, pixel p ) -- places pixel at (x,y)
 *
 */
#include <cstdlib>
#include <cstring>
#include <ctype.h>
typedef unsigned char BYTE;
struct pixel_struct { int r, g, b; };
struct image_struct {
	int width, height, max, bpp, raster;
	BYTE * raw;
	BYTE ** pixel;    
};
typedef struct image_struct image;
typedef struct pixel_struct pixel;

pixel get_pixel( image z, int x, int y )
{
	BYTE * foo;
	pixel p;
	foo = z.pixel[y]+(x*3*z.bpp);
	if( z.bpp==1 ) { p.r = foo[0]; p.g = foo[1]; p.b = foo[2]; }
	else {
		p.r = (foo[0]<<8) + foo[1];
		p.g = (foo[2]<<8) + foo[3];
		p.b = (foo[4]<<8) + foo[5];
	}
	return p;
}

void put_pixel( image z, int x, int y, pixel p )
{
	BYTE * foo;
	foo = z.pixel[y]+(x*3*z.bpp);
	if( z.bpp==1 ) { foo[0] = p.r; foo[1] = p.g; foo[2] = p.b; }
	else {
		foo[0] = (p.r>>8) & 0xff;
		foo[1] = p.r & 0xff;
		foo[2] = (p.g>>8) & 0xff;
		foo[3] = p.g & 0xff;
		foo[4] = (p.b>>8) & 0xff;
		foo[5] = p.b & 0xff;
	}
}

int alloc_ppm( image * y ) {
	int i=0;
	y->raw    = (BYTE *) malloc( y->height * y->raster );
	y->pixel  = (BYTE **) malloc( y->height * sizeof(BYTE*) );
	if( y->raw!=NULL && y->pixel!=NULL ) 
		for( i=0; i<y->height; i++ ) 
			y->pixel[i] = y->raw + i*y->raster;
	else fprintf( stderr, "Failed to allocate memory!\n" );
	return i;
}

void dealloc_ppm( image x ) {
	free( x.pixel );
	free( x.raw );
}

int copy_ppm( image x, image * y ) {
	int i=0;
	*y = x;
	if( alloc_ppm(y) ) 
		for( i=0; i<y->height*y->raster; i++ ) 
			y->raw[i] = x.raw[i];
	return i;
}

int get_int( void ) { /* scan int from ppm, passing over comments */
	int k;
	unsigned char ch;
	do {
		while( isspace(ch=getchar()) ) ;
		if( ch=='#' ) while( (ch=getchar())!='\n' ) ;
	} while( !isdigit(ch) );
	for( k=(ch-'0'); isdigit(ch=getchar()); k = 10*k+(ch-'0') );
	return k;
}

int get_ppm( image * x) {
	int flag=0;
	char magic[3]="ZZ";

	fread(magic,2,1,stdin);	
	if( strcmp(magic,"P6")==0 ) {
		x->width=get_int(); x->height=get_int(); x->max=get_int();
		x->bpp = (x->max<256) ? 1 : 2;
		x->raster = x->width * x->bpp * 3;
		if( flag=alloc_ppm(x) ) 
			fread( x->raw, x->height, x->raster, stdin );
	}
	else fprintf( stderr, "This is not a valid binary PPM (magic # is not P6)\n" );
	return flag;
}

void put_ppm( image x )
{
	printf( "P6\n%d %d %d\n", x.width, x.height, x.max ) ;
	fwrite( x.raw, x.height, x.raster, stdout );
}
