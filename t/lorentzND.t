#!perl

use PDL;
use PDL::MatrixOps;
use PDL::Transform;

use Math::Trig qw[ pi ];

use Test::More tests => 8;

BEGIN {
  use_ok('PDL::FuncND');
}

use strict;
use warnings;

our $N = 10;

sub _lorentz1D {

    my ( $x, $c, $s ) = @_;

    return $s**2 / ( ($x - $c)**2 + $s**2);

}

sub _lorentz2D_symmetric {

    my ( $x, $y, $c, $s ) = @_;

    return $s**2 / ( ( $x - $c->[0] )**2 + ( $y - $c->[1] )**2 + $s**2);

}



# Perl
sub _lorentz2D {

    my ( $x, $y, $c, $s ) = @_;

    my ($x0, $y0) = @$c;
    my ( $s00, $s01, $s10, $s11 ) = $s->list;


    return 1/(($s00*($y-$y0)/($s00*$s11-$s01*$s10)
	       - $s10*($x-$x0)/($s00*$s11-$s01*$s10))
	      * ($y-$y0)
	      + ($x-$x0)
	      * ($s11*($x-$x0)/($s00*$s11-$s01*$s10)
		 - $s01*($y-$y0)/($s00*$s11-$s01*$s10))+1);

}

# test 1D
{
    my $c = 2.2;
    my $x = sequence(10);
    my $s = 3;


    my $f = _lorentz1D( $x, $c, $s );
    my $g = lorentzND( $x, { scale => $s, center => $c } );

    ok( all ( approx($g, $f) ), "1D" );
}

# test 1D, transformed
{
    my $c = 2.2;
    my $x = 1.3 * ( sequence(10) + 1.9 );
    my $s = 3;

    my $t = t_linear( { scale => 1.3, offset => 1.9, dims => 1 } );

    my $f = _lorentz1D( $x, $c, $s );
    my $g = lorentzND( $x, { scale => $s, center => $c, transform => $t } );

    ok( all ( approx($g, $f) ), "1D, transformed" );
}


# test 2D, symmetric
{

    my @c = ( 2.2, 9.3);
    my $x = xvals(10,10);
    my $y = yvals(10,10);

    my $s = 3;

    my $f = _lorentz2D_symmetric( $x, $y, \@c, $s );

    my $g = lorentzND( $x, { scale => $s, center => \@c } );

    ok( all ( approx($g, $f) ), q[2D, symmetric] );


    $g = lorentzND( $x, { scale => [$s,$s], center => \@c } );

    ok( all ( approx($g, $f) ), q[2D, symmetric, array] );


}

# test 2D, symmetric, transformed
{
    my @o = ( 8.3, 7.7 );
    my @m = ( 1.1, 3.9 );

    my @c = ( 2.2, 9.3);
    my $x = $m[0] * ( xvals(10,10) + $o[0] );
    my $y = $m[1] * ( yvals(10,10) + $o[1] );

    my $t = t_linear( { scale => pdl( @m ), offset => pdl( @o ) } );

    my $s = 3;

    my $f = _lorentz2D_symmetric( $x, $y, \@c, $s );

    my $g = lorentzND( $x, { scale => $s, center => \@c, transform => $t } );

    ok( all ( approx($g, $f) ), "2D symmetric, transformed" );
}


# test 2D, general
{

    my @c = ( 2.2, 9.3);
    my $x = xvals(10,10);
    my $y = yvals(10,10);

    my $s = pdl( [1,2], [3,8] );

    my $f = _lorentz2D( $x, $y, \@c, $s );

    my $g = lorentzND( $x, { scale => $s, center => \@c } );

    ok( all ( approx($g, $f) ), q[2D] );


}


# test 2D, general, transformed
{
    my @o = ( 8.3, 7.7 );
    my @m = ( 1.1, 3.9 );

    my @c = ( 2.2, 9.3);
    my $x = $m[0] * ( xvals(10,10) + $o[0] );
    my $y = $m[1] * ( yvals(10,10) + $o[1] );

    my $s = pdl( [1,2], [3,8] );

    my $t = t_linear( { scale => pdl( @m ), offset => pdl( @o ) } );

    my $f = _lorentz2D( $x, $y, \@c, $s );

    my $g = lorentzND( $x, { scale => $s, center => \@c, transform => $t } );

    ok( all ( approx($g, $f) ), "2D general, transformed" );
}

