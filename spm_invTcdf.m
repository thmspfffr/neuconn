function x = spm_invTcdf(F,v)
2	
% Inverse Cumulative Distribution Function (CDF) of Students t distribution
3	
% FORMAT x = spm_invTcdf(F,v)
4	
%
5	
% F  - CDF (lower tail p-value)
6	
% v - degrees of freedom (v>0, non-integer d.f. accepted)
7	
% x - T-variate (Student's t has range (-Inf,Inf))
8	
%__________________________________________________________________________
9	
%
10	
% spm_invTcdf implements the inverse Cumulative Distribution of Students
11	
% t-distributions.
12	
%
13	
% Definition:
14	
%--------------------------------------------------------------------------
15	
% The Student's t-distribution with v degrees of freedom is defined for
16	
% positive integer v and real x. The Cumulative Distribution
17	
% Function (CDF) F(x) is the probability that a realisation of a
18	
% t-distributed random variable X has value less than x. F(x)=Pr{X<x}:
19	
% (See Evans et al., Ch37)
20	
%
21	
% This implementation is not restricted to whole (positive integer) df
22	
% v, rather it will compute for any df v>0.
23	
%
24	
% Variate relationships: (Evans et al., Ch37 & 7)
25	
%--------------------------------------------------------------------------
26	
% The Student's t distribution with 1 degree of freedom is the Standard
27	
% Cauchy distribution, which has a simple closed form CDF.
28	
%
29	
% For X a t-variate with v degrees of freedom, v/(v+X^2) has
30	
% distribution related to a Beta random variable with shape parameters
31	
% w/2 & 1/2, as described below.
32	
%
33	
% Algorithm:
34	
%--------------------------------------------------------------------------
35	
% Using the routine spm_invBcdf for the Beta distribution, with
36	
% appropriate parameters:  The CDF of the Student's t-distribution with
37	
% v degrees of freedom is related to the incomplete beta function by:
38	
%       Pr(|X|<x) = betainc(v/(v+x^2),v/2,1/2)
39	
% so
40	
%              {     betainc(v/(v+x^2),v/2,1/2) / 2      for x<0
41	
%       F(x) = |   0.5                                   for x=0
42	
%              { 1 - betainc(v/(v+x^2),v/2,1/2) / 2      for x>0
43	
%
44	
% See Abramowitz & Stegun, 26.5.27 & 26.7.1; Press et al., Sec6.4 for
45	
% definitions of the incomplete beta function. The relationship is
46	
% easily verified by substituting for v/(v+x^2) in the integral of the
47	
% incomplete beta function.
48	
%
49	
% References:
50	
%--------------------------------------------------------------------------
51	
% Evans M, Hastings N, Peacock B (1993)
52	
%       "Statistical Distributions"
53	
%        2nd Ed. Wiley, New York
54	
%
55	
% Abramowitz M, Stegun IA, (1964)
56	
%       "Handbook of Mathematical Functions"
57	
%        US Government Printing Office
58	
%
59	
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
60	
%       "Numerical Recipes in C"
61	
%        Cambridge
62	
%
63	
%__________________________________________________________________________
64	
% Copyright (C) 1993-2011 Wellcome Trust Centre for Neuroimaging
65	
66	
% Andrew Holmes
67	
% $Id: spm_invTcdf.m 4182 2011-02-01 12:29:09Z guillaume $
68	
69	
70	
%-Format arguments, note & check sizes
71	
%--------------------------------------------------------------------------
72	
if nargin<2, error('Insufficient arguments'), end
73	
74	
ad = [ndims(F);ndims(v)];
75	
rd = max(ad);
76	
as = [  [size(F),ones(1,rd-ad(1))];...
77	
    [size(v),ones(1,rd-ad(2))]     ];
78	
rs = max(as);
79	
xa = prod(as,2)>1;
80	
if all(xa) && any(diff(as(xa,:)))
81	
    error('non-scalar args must match in size');
82	
end
83	
84	
85	
%-Computation
86	
%--------------------------------------------------------------------------
87	
%-Initialise result to zeros
88	
x = zeros(rs);
89	
90	
%-Only defined for F in [0,1] & strictly positive v.
91	
% Return NaN if undefined.
92	
md = ( F>=0  &  F<=1  &  v>0 );
93	
if any(~md(:))
94	
    x(~md) = NaN;
95	
    warning('Returning NaN for out of range arguments');
96	
end
97	
98	
%-Special case: x is 0 when F=0.5, -Inf when F=0, +Inf when F=1
99	
x(md & F==0) = -Inf;
100	
x(md & F==1) = +Inf;
101	
102	
%-Special case: Standard Cauchy distribution when v=1
103	
ml = ( md  &  v==1 ); if xa(1), mlF=ml; else mlF=1; end
104	
x(ml) = tan(pi*(F(mlF)-0.5));
105	
106	
%-Compute where defined & not special cases
107	
Q  = find( md  &  F>0  &  F~=0.5  &  F<1  &  v~=1 );
108	
if isempty(Q), return, end
109	
if xa(1), QF=Q; else QF=1; end
110	
if xa(2), Qv=Q; else Qv=1; end
111	
112	
%-Compute
113	
xQPos = F(QF)>0.5;
114	
bQ    = spm_invBcdf(2*(xQPos -(xQPos*2-1).*F(QF)),v(Qv)/2,1/2);
115	
x(Q)  = (xQPos*2-1) .* sqrt(v(Qv)./bQ -v(Qv));