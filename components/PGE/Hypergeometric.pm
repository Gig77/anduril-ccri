package Hypergeometric;

# Hypergeometric test implementation
# Note: in case numbers are too big, pbinom from Math::CDF is returned
#
# Copyright (C) 2007 Roland Barriot
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#		
# email: 
#	roland.barriot@gmail.com
# mail: 
# 	Roland Barriot
# 	ESAT-SCD
# 	Katholieke Universiteit Leuven
# 	Kasteelpark Arenberg, 10
# 	B-3001 Leuven
# 	Belgium

use Math::CDF;

# p-value(s,q,t,c) = probability to observe at leat c common elements between to random samples of size q and t drawn on the same population of size s.

sub new {
    my $class = shift;
    my $self = {};
    bless ($self, $class);
    return $self;
}


sub compute {
    my $self = shift;
    my $populationSize = shift;
    my $querySet = shift;
    my $targetSet = shift;
    my $commonElements = shift;
    return hypergeometric($populationSize,$querySet,$targetSet,$commonElements);
}

sub hypergeometric {
    my $s = shift;
    my $q = shift;
    my $t = shift;
    my $c = shift;
    my $res = 0;
    my $min=$q;
    if ($t < $q) { $min=$t; }
    my $cnp_s_q = cnp($s,$q);

    if ($cnp_s_q eq 'inf') {
	return Math::CDF::pbinom($q - $c,               $q,      1 - $t/$s);
    }
    for (my $i=$c ; $i <= $min ; $i++) {
	$res += (cnp($t,$i) * cnp($s-$t, $q-$i) / $cnp_s_q);
    }
    return $res;
}

# / n \
# |   |
# \ p /
sub cnp {
    my $n = shift;
    my $p = shift;
    if    ($p > $n)  { return 0; }
    elsif ($p == 0)  { return 1; }
    elsif ($p == 1)  { return $n; }
    elsif ($p == $n) { return 1; }
    else {
	if ($p > $n/2) { $p = $n-$p; }
	my $res=$n; # p==1 -> Cnp=n
	for (my $i=2 ; $i <= $p ; $i++) {
	    $res = ($res * ($n-$i+1))/$i;
	}
	return $res;
    }    
}


1;
