# --8<--8<--8<--8<--
#
# Copyright (C) 2010 Smithsonian Astrophysical Observatory
#
# This file is part of PDL::FuncND
#
# PDL::FuncND is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -->8-->8-->8-->8--

package PDL::FuncND;

use strict;
use warnings;

use Math::Trig qw[ pi ];

use PDL;
use PDL::Math;
use PDL::Transform qw[];
use PDL::Options qw[ iparse ];

use Math::Trig qw[ pi ];

use base 'PDL::Exporter';

our @EXPORT_OK = qw[
		       mahalanobis
		       cauchyND
		       gaussND
		       lorentzND
		  ];

our %EXPORT_TAGS = (Func     => [@EXPORT_OK],
		   );

# bizarre spurious errors are being thrown by this policy
## no critic (ProhibitAccessOfPrivateData)

our $VERSION = '0.06';

=head1 NAME

PDL::FuncND - N dimensional version of functions


=head1 SYNOPSIS

    use PDL::FuncND;


=head1 DESCRIPTION

This module provides multi-dimensional implementations of common
functions.

=head1 FUNCTIONS

=cut


# keep option handling uniform.  barf if there's a problem
sub _handle_options {


    # remove known arguments from @_ so can use new_from_specification
    my $self = shift;
    my $opt  = shift;

    my ( $vectors, $output, $ndims, $center, $scale );

    if ( $opt->{vectors} )
    {
	barf( "first argument must be a piddle if vectors option is set\n" )
	  unless eval { $self->isa( 'PDL' ) };

	$output = shift;
	$vectors = $self;

	barf( "wrong dimensions for input vectors; expected 2 , got ",
	      $vectors->ndims, "\n" )
	  unless $vectors->ndims == 2;

	( $ndims, my $m ) = $vectors->dims;

	$vectors = $opt->{transform}->apply( $vectors )
	  if defined $opt->{transform};

	if ( defined $opt->{center} )
	{
	    barf( "cannot use center = 'auto' if vectors is set\n" )
	      if $opt->{center} eq 'auto';
	    $center = PDL::Core::topdl( $opt->{center} );
	}
    }

    else
    {
	$output = @_ ? $self->new_from_specification( @_ )
		     : $self->new_or_inplace;

	$ndims = $output->ndims;

	$vectors = $output->ndcoords->reshape( $ndims, $output->nelem );

	$vectors->inplace->apply( $opt->{transform} )
	  if defined $opt->{transform};

	if ( defined $opt->{center} )
	{
	    if ( $opt->{center} eq 'auto' )
	    {
		$center = ( pdl([$output->dims]) - 1 ) / 2;
		$center->inplace->apply( $opt->{transform} )
		  if defined $opt->{transform};
	    }
	    else
	    {
		$center =  PDL::Core::topdl( $opt->{center} )
	    }
	}
    }

    # for 1D output $center may be a 0D PDL; this causes grief;
    # promote it to a 1D
    $center = $center->dummy(0) if defined $center && $center->ndims == 0;

    barf( "center vector has wrong dimensions\n" )
      if defined $center && ($center->ndims != 1
			     || ($center->dims)[0] != $ndims);

    # handle scale
    # scalar -> symmetric, independent
    if ( ! ref $opt->{scale} )
    {
	$scale = identity( $ndims ) * $opt->{scale} ** 2;
    }

    # 1D piddle of length N
    elsif ( 'ARRAY' eq ref $opt->{scale}
	    && @{$opt->{scale}} == $ndims
	  )
    {
	$scale = identity( $ndims ) * pdl( $opt->{scale} ) ** 2;
    }

    # 1D piddle of length N
    elsif ( eval { $opt->{scale}->isa('PDL') }
	    && $opt->{scale}->ndims == 1
	    && ($opt->{scale}->dims)[0] == $ndims
	  )
    {
	$scale = identity( $ndims ) * $opt->{scale} ** 2;
    }

    # full  matrix N^N piddle
    elsif ( eval { $opt->{scale}->isa('PDL') }
	    && $opt->{scale}->ndims == 2
	    && all( pdl( $opt->{scale}->dims ) == pdl( $ndims, $ndims ) )
	  )
    {
	$scale = $opt->{scale};
    }

    else
    {
	barf( "scale argument is not a scalar, an array of length $ndims,",
	      " or a piddle of shape ($ndims) or shape ($ndims,$ndims)\n" );
    }

    # apply a rotation to the scale matrix
    if ( defined $opt->{theta} ) {
	barf( "theta may only be used for 2D PDFs\n" )
	  if $ndims != 2;

	my $R = pdl( [ cos($opt->{theta}), -sin($opt->{theta} ) ],
		     [ sin($opt->{theta}),  cos($opt->{theta} ) ] );
	$scale = $R->transpose x $scale x $R;
    }



    return ( vectors => $vectors,
	     output  => $output,
	     ndims   => $ndims,
	     center  => $center,
	     ndims   => $ndims,
	     scale   => $scale
	   );



}


=pod

=head2 cauchyND

=for ref

Evaluate the multi-variate Cauchy function on an N-dimensional grid or
at a set of locations.


=for usage

  $a = cauchyND( [OPTIONAL TYPE], $nx, $ny, ..., \%options );
  $b = cauchyND( $a, \%options );
  cauchyND( inplace $a, \%options );
  $a->inplace->cauchyND( \%options );


B<cauchyND> can evaluate the function either on a grid or at discrete
locations:

=over

=item * evaluation on a grid

Either specify the output piddle dimensions explicitly,

  $f = cauchyND( [ OPTIONAL TYPE], $nx, $ny, ..., \%options );

or specify a template piddle I<without> specifying the C<vectors>
option:

  $f = cauchyND( $piddle, \%options );

By default B<cauchyND> will evaluate the function at the I<indices> of
the points in the input piddle.  These may be mapped to other values
by specifying a transform with the C<transform> option.  B<cauchyND>
is inplace aware, and will use B<$piddle> as the output piddle if its
inplace flag is set.

  cauchyND( inplace $f, \%options );
  $f->inplace->cauchyND( \%options );

=item * evaluation at a set of locations

The input piddle should represent a set of vectors and should have a
shape of (N,m), where C<m> is the number of vectors in the set. The
C<vectors> option must also be set:

  $piddle = pdl( [2,1], [3,1], [4,2]  );
  $f = cauchyND( $piddle, { vectors => 1 } );

The vectors may be transformed before use via the C<transform> option.

=back

The following options are available:

=over

=item C<center> | C<centre>

The center of the distribution.  If not specified it defaults to the
origin.

This may take one of the following values:

=over

=item * A vector of shape (N).

The location of the center. This may be either a Perl arrayref or a
one dimensional piddle.  If the input coordinates are transformed,
this is in the I<transformed> space.

=item * the string C<auto>

If the PDF is calculated on a grid, this will center the distribution on
the grid. It is an error to use this for explicit locations.

=back

=item C<transform>

A PDL::Transform object to be applied to the input coordinates.

=item C<scale>

The scale. If the input coordinates are transformed
via the C<transform> option, the units of scale are those in the
I<transformed> space.  This may be specified as:

=over

=item * a scalar (Perl or piddle)

This results in a symmetric distribution with the given scale along each
coordinate.

=item * a vector of shape (N) (piddle or Perl arrayref)

This results in a distribution with the specified scales for each
coordinate.

=item * a matrix (piddle of shape (N,N))

This should be a positive-definite matrix containing squared
scales.

=back

=item C<theta> (Perl scalar)

B<Only for 2D!> Applies a rotation (clockwise, e.g. +X
rotates towards -Y) by the specified angle (specified in radians).

=item C<log> (Boolean)

If true, return the logarithm of the function. Defaults to false.

=back

=cut

# from http://en.wikipedia.org/wiki/Cauchy_distribution#Multivariate_Cauchy_distribution

#                                    1 + k
#                              Gamma(-----)
#                                      2
#     ------------------------------------------------------------
#                                                          1 + k
#                                                          -----
#           1     k/2    1/2              T   -1             2
#     Gamma(-)  pi    |S|    (1 + (x - mu)   S    (x - mu))
#           2
#


sub _gamma { exp( (lgamma(@_))[0]) }

sub cauchyND {

    # handle being called as a method or a function
    my $self = eval { ref $_[0] && $_[0]->isa( 'PDL' ) } ? shift @_ : 'PDL';

    # handle options.
    my $opt = 'HASH' eq ref $_[-1] ? pop( @_ ) : {};
    my %opt = iparse( { center    => undef,
			scale     => 1,
			transform => undef,
			vectors   => 0,
			log       => 0,
			theta     => undef,
		      }, $opt );

    my %par = _handle_options( $self, \%opt, @_ );

    my $d2 = mahalanobis( $par{vectors}, $par{scale},
			    { squared => 1,
			      ( defined $par{center} ? (center => $par{center}) : () )
			    });

    my $k = $par{ndims};

    my $pdf =
      _gamma((1+$k)/2)
      /  (   _gamma(1/2)
	   * pi**($k/2)
	   * sqrt(determinant( $par{scale} ))
	   * (1 + $d2) ** ( (1+$k)/2 )
	 );

    my $retval = $opt{log} ? log($pdf) : $pdf;

    my $output = $par{output};

    if ( $opt{vectors} )
    {
	$output = $retval;
    }

    else
    {
	$output .= $retval->reshape( $output->dims );
    }

    return $output;
}

*PDL::cauchyND = \&cauchyND;


=pod

=head2 gaussND

=for ref

Evaluate the sampled multi-dimensional Gaussian PDF on an
N-dimensional grid or at a set of locations.

=for usage

  $f = gaussND( [OPTIONAL TYPE], $nx, $ny, ..., \%options );
  $f = gaussND( $piddle, \%options );
  gaussND( inplace $piddle, \%options );
  $a->inplace->gaussND( \%options );


B<gaussND> can evaluate the function either on a grid or at discrete
locations:

=over

=item * evaluation on a grid

Either specify the output piddle dimensions explicitly,

  $f = gaussND( [ OPTIONAL TYPE], $nx, $ny, ..., \%options );

or specify a template piddle I<without> specifying the C<vectors>
option:

  $f = gaussND( $piddle, \%options );

By default B<gaussND> will evaluate the function at the I<indices> of
the points in the input piddle.  These may be mapped to other values
by specifying a transform with the C<transform> option.  B<gaussND> is
inplace aware, and will use B<$piddle> as the output piddle if its
inplace flag is set.

  gaussND( inplace $f, \%options );
  $f->inplace->gaussND( \%options );

=item * evaluation at a set of locations

The input piddle should represent a set of vectors and should have a
shape of (N,m), where C<m> is the number of vectors in the set. The
C<vectors> option must also be set:

  $piddle = pdl( [2,1], [3,1], [4,2]  );
  $f = gaussND( $piddle, { vectors => 1 } );

The vectors may be transformed before use via the C<transform> option.

=back

The following options are available:

=over

=item C<center> | C<centre>

The center of the distribution.  If not specified it defaults to the
origin.

This may take one of the following values:

=over

=item * A vector of shape (N).

The location of the center. This may be either a Perl arrayref or a
one dimensional piddle.  If the input coordinates are transformed,
this is in the I<transformed> space.

=item * the string C<auto>

If the PDF is calculated on a grid, this will center the distribution on
the grid. It is an error to use this for explicit locations.

=back

=item C<transform>

A PDL::Transform object to be applied to the input coordinates.

=item C<scale>

The scale. If the input coordinates are transformed
via the C<transform> option, the units of scale are those in the
I<transformed> space.  This may be specified as:

=over

=item * a scalar (Perl or piddle)

This results in a symmetric distribution with the given scale along each
coordinate.

=item * a vector of shape (N) (piddle or Perl arrayref)

This results in a distribution with the specified scales for each
coordinate.

=item * the full covariance matrix (piddle of shape (N,N))

This results in a distribution with correlated scales. At present this
matrix is not verified to be a legitimate covariance matrix.

=back

=item C<theta> (Perl scalar)

B<Only for 2D!> Applies a rotation (clockwise, e.g. +X
rotates towards -Y) by the specified angle (specified in radians).

=item C<log> (Boolean)

If true, return the logarithm of the function. Defaults to false.

=back

=cut


sub gaussND {

    # handle being called as a method or a function
    my $self = eval { ref $_[0] && $_[0]->isa( 'PDL' ) } ? shift @_ : 'PDL';

    # handle options.
    my $opt = 'HASH' eq ref $_[-1] ? pop( @_ ) : {};
    my %opt = iparse( { center    => undef,
			scale     => 1,
			transform => undef,
			vectors   => 0,
			log       => 0,
			theta     => undef,
		      }, $opt );

    my %par = _handle_options( $self, \%opt, @_ );

    my $d2 = mahalanobis( $par{vectors}, $par{scale},
			{ squared => 1,
			  ( defined $par{center} ? (center => $par{center}) : () )
			});
    my $log_pdf = -(   $par{ndims} * log(2 * pi)
		     + log( determinant( $par{scale} ))
		     + $d2
		   )/2;

    my $retval = $opt{log} ? $log_pdf : exp( $log_pdf );

    my $output = $par{output};

    if ( $opt{vectors} )
    {
	$output = $retval;
    }

    else
    {
	$output .= $retval->reshape( $output->dims );
    }

    return $output;
}

*PDL::gaussND = \&gaussND;


=pod

=head2 lorentzND

=for ref

Evaluate the multi-dimensional Lorentz function on an
N-dimensional grid or at a set of locations.

=for usage

  $f = lorentzND( [OPTIONAL TYPE], $nx, $ny, ..., \%options );
  $f = lorentzND( $piddle, \%options );
  lorentzND( inplace $piddle, \%options );
  $a->inplace->lorentzND( \%options );


The Lorentz function is usually defined in one dimension as.

                       2
                      g
  f(x; x0, g) = --------------
                       2    2
                (x - x0)  + g


where I<g> is the half-width at half-max (HWHM).  The two dimensional
symmetric analogue (sometimes called the "radial Lorentz
function") is

                                    2
                                   g
  f(x, y; x0, y0, g) = --------------------------
                               2           2    2
                       (x - x0)  + (y - y0)  + g


One can extend this to an asymmetric form by writing it as

                            1
  f(x; u, S) = ---------------------------
                      T    -1
               (x - u)  . S  . (x - u) + 1

where I<x> is now a vector, I<u> is the expectation value of the
distribution, and I<S> is a matrix describing the N-dimensional scale
of the distribution akin to (but not the same as!) a covariance matrix.

For example, a symmetric 2D Lorentz with HWHM of I<g> has

       [  2     ]
       [ g   0  ]
  S =  [        ]
       [      2 ]
       [ 0   g  ]

while an elliptical distribution elongated twice as much along the
I<X> axis as the I<Y> axis would be:

       [     2      ]
       [ (2*g)   0  ]
  S =  [            ]
       [          2 ]
       [ 0       g  ]


B<lorentzND> extends the Lorentz function to N dimensions by treating
I<x> and I<u> as vectors of length I<N>, and I<S> as an I<NxN> matrix.

It can evaluate the function either on a grid or at discrete
locations:

=over

=item * evaluation on a grid

Either specify the output piddle dimensions explicitly,

  $f = lorentzND( [ OPTIONAL TYPE], $nx, $ny, ..., \%options );

or specify a template piddle I<without> specifying the C<vectors>
option:

  $f = lorentzND( $piddle, \%options );

By default B<lorentzND> will evaluate the function at the I<indices> of
the points in the input piddle.  These may be mapped to other values
by specifying a transform with the C<transform> option.  B<lorentzND> is
inplace aware, and will use B<$piddle> as the output piddle if its
inplace flag is set.

  lorentzND( inplace $f, \%options );
  $f->inplace->lorentzND( \%options );

=item * evaluation at a set of locations

The input piddle should represent a set of vectors and should have a
shape of (N,m), where C<m> is the number of vectors in the set. The
C<vectors> option must also be set:

  $piddle = pdl( [2,1], [3,1], [4,2]  );
  $f = lorentzND( $piddle, { vectors => 1 } );

The vectors may be transformed before use via the C<transform> option.

=back

The following options are available:

=over

=item C<center> | C<centre>

The center of the distribution.  If not specified it defaults to the
origin.

This may take one of the following values:

=over

=item * A vector of shape (N).

The location of the center. This may be either a Perl arrayref or a
one dimensional piddle.  If the input coordinates are transformed,
this is in the I<transformed> space.

=item * the string C<auto>

If the PDF is calculated on a grid, this will center the distribution on
the grid. It is an error to use this for explicit locations.

=back

=item C<transform>

A PDL::Transform object to be applied to the input coordinates.

=item C<scale>

The scale. If the input coordinates are transformed
via the C<transform> option, the units of scale are those in the
I<transformed> space.  This may be specified as:

=over

=item * a scalar (Perl or piddle)

This results in a symmetric distribution with the given scale along each
coordinate.

=item * a vector of shape (N) (piddle or Perl arrayref)

This results in a distribution with the specified scales for each
coordinate.

=item * a matrix (piddle of shape (N,N))

This should be a positive-definite matrix containing squared
scales.

=back

=item C<theta> (Perl scalar)

B<Only for 2D!> Applies a rotation (clockwise, e.g. +X
rotates towards -Y) by the specified angle (specified in radians).

=back

=cut

sub lorentzND {

    # handle being called as a method or a function
    my $self = eval { ref $_[0] && $_[0]->isa( 'PDL' ) } ? shift @_ : 'PDL';

    # handle options.
    my $opt = 'HASH' eq ref $_[-1] ? pop( @_ ) : {};
    my %opt = iparse( { center    => undef,
			scale     => 1,
			transform => undef,
			vectors   => 0,
			theta     => undef,
		      }, $opt );

    my %par = _handle_options( $self, \%opt, @_ );

    my $d2 = mahalanobis( $par{vectors}, $par{scale},
			{ squared => 1,
			  ( defined $par{center} ? (center => $par{center}) : () )
			});
    my $f = 1/ ( 1 + $d2 );

    my $retval = $f;

    my $output = $par{output};

    if ( $opt{vectors} )
    {
	$output = $retval;
    }

    else
    {
	$output .= $retval->reshape( $output->dims );
    }

    return $output;
}

*PDL::lorentzND = \&lorentzND;



=pod

=head2 mahalanobis

=for ref

Calculate the Mahalanobis distance for one or more vectors

=for sig

  Signature: ( x(n,m), s(n,n), [o]d(m), \%options )

=for usage

  $d = mahalanobis( $v, $S, \%options );
  mahalanobis( $v, $S, $d, \%options );

The Mahalanobis distance of a multivariate vector (v) from a location
(u) with a covariance matrix (S) is defined as

  dm(x,u) = sqrt( (v-u)T S^-1 (v-u) )

The input piddle representing the vectors (C<$v>) must have shape (N,m),
where C<N> is the dimension of the vector space and C<m> is the number
of vectors.

The input covariance matrix (C<$S>) must have shape (N,N).  It is I<not>
checked for validity.

The available options are:

=over

=item C<center> | C<centre>

The vector from which the distance is to be calculated.  It must have shape (N).
It defaults to the origin.

=item C<inverted>

If true, the input matrix is the inverse of the covariance matrix.
Defaults to false.

=item C<squared>

if true, the returned values are the distances squared.

=back


=cut

sub mahalanobis {

    # handle options.
    my $opt = 'HASH' eq ref $_[-1] ? pop( @_ ) : {};
    my %opt = PDL::Options::iparse( { center => undef,
				      inverted => 0,
				      squared  => 0,
				    }, $opt );

    my ( $x, $scale, $out ) = @_;

    # invert the matrix if it hasn't already been inverted
    $scale = $scale->inv
	unless $opt{inverted};

    my $xc;
    if ( defined $opt{center} )
    {
	my $c = PDL::Core::topdl( $opt{center} );
	$xc = $x - $c;
    }
    else
    {
	$xc = $x;
    }

    $out = zeroes(double, ($x->dims)[-1] )
      unless defined $out;

    inner2( $xc, $scale, $xc, $out );

    $out->inplace->sqrt unless $opt{squared};

    return $out;
}

*PDL::mahalanobis = \&mahalanobis;


__END__



=head1 SEE ALSO

PDL::Func.

=head1 BUGS

Please report bugs to https://rt.cpan.org/Public/Dist/Display.html?Name=PDL-FuncND.

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010-2012 The Smithsonian Astrophysical Observatory

PDL::FuncND is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 AUTHOR

Diab Jerius  E<lt>djerius@cpan.orgE<gt>
