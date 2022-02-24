#include <stdio.h>
#include "ppm.c"

int main( void )   /* usage:  foo.exe < in.ppm > out.ppm */
{
	image in, out;
	pixel p;
	int x,y;

	if( get_ppm( &in ) && copy_ppm( in, &out ) ) {

		for( x=0; x<in.width; x++ ) 
			for( y=0; y<in.height; y++ ) {
				p = get_pixel( in, x, y );			
				p.r = 0;
				p.b = 0;
				put_pixel( out, x, y, p );
			}

		put_ppm( out );
		dealloc_ppm( in );
		dealloc_ppm( out );
	}
	else
	{
		printf("hello\n");
	}
	return 0;
}
