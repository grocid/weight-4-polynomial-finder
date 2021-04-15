# Weight-4 Polynomial Multiple Finder

## Finding multiples of polynomials over 𝔽<sub>2</sub> 

The problem of finding a low-weight multiple of a polynomial is of great importance in the field of cryptography. Some applications are distinguishing and correlations attacks on stream ciphers. Another area is finite-field arithmetic.

We define problem as follows. Given a polynomial *P*(*x*), find another polynomial *u*(*x*) such that the resulting product has low weight (or consists of few terms).

It is not known exactly how hard the problem over 𝔽<sub>2</sub> is, but it is believed to be hard. One simple reduction is to represent the polynomial as a linear code and finding a minimum-weight codeword using information-set decoding. Other methods such as time--memory trade-offs are applicable as well as *k*-sum algorithms. Also memory-less constructions are possible. Such an algorithm is possible by defining two functions *f* and *g*, that act in any input as a random oracle to ℤ<sub>m</sub>, where *m* is close to the degree of *K*(*x*). Then, a cycle-finding algorithm can be run for the recursion *y* = *x*<sup>*f*(*y*)</sup> + *x*<sup>*g*(y)</sup>. When the two runners collide, there are two values *y* and *y*' such that *x*<sup>*f*(*y*)</sup> + *x*<sup>*g*(*y*)</sup> = *x*<sup>*f*(*y*')</sup> + *x*<sup>*g*(*y*')</sup> which gives the exponents for the multiple *K*(*x*) if *y* ≠ *y*'.

For the case *w* = 4, a 4-sum algorithm (generalized-birthday algorithm) is used. By extending the range from 2<sup>*d*/3</sup> to 2<sup>*d*/3 + _α_</sup>, the probability increases considerably. To see this, assume that *K*(*x*) = *x*<sup>i</sup> + *x*<sup>j</sup> + *x*<sup>k</sup> + *x*<sup>l</sup> has degree *m* = 2<sup>*d*/3</sup>. By the extension of the range, there are *N* = 2<sup>*d*/3 + α</sup> - 2<sup>*d*/3</sup> = 2<sup>*d*/3</sup> ( 2<sup>α</sup> - 1) shifts of the form *x*<sup>i</sup> *K*(*x*) in this set. Using a mask (defined as the function φ) of weight *d*/3, the probability that φ (*x*<sup>i</sup> + *x*<sup>j</sup>) = 0 is  2<sup>-*d*/3</sup>, and φ  (*x*<sup>k</sup> + *x*<sup>l</sup>) = 0 follows by definition. Since this is repeated for all *N* shifts, the probability that one passes through the filter is overwhelmingly large.

## Implementation
The implementation is a threaded and scalable to up arbitrarily many cores. It is also constructed in a way such that mutexes are not needed.

## Results


Origin | Polynomial | Multiple
--- | --- | ---
Bluetooth Core E0 | *x*<sup>97</sup> + *x*<sup>92</sup> + *x*<sup>89</sup> + *x*<sup>80</sup> + *x*<sup>76</sup> + *x*<sup>73</sup> + *x*<sup>69</sup> + *x*<sup>68</sup> + *x*<sup>64</sup> + *x*<sup>61</sup> + *x*<sup>59</sup> + *x*<sup>58</sup> + *x*<sup>57</sup> + *x*<sup>56</sup> + *x*<sup>55</sup> + *x*<sup>54</sup> + *x*<sup>53</sup> + *x*<sup>51</sup> + *x*<sup>49</sup> + *x*<sup>44</sup> + *x*<sup>43</sup> + *x*<sup>41</sup> + *x*<sup>39</sup> + *x*<sup>36</sup> + *x*<sup>35</sup> + *x*<sup>34</sup> + *x*<sup>32</sup> + *x*<sup>31</sup> + *x*<sup>30</sup> + *x*<sup>21</sup> + *x*<sup>19</sup> + *x*<sup>15</sup> + *x*<sup>13</sup> + *x*<sup>9</sup> + *x*<sup>5</sup> + *x*<sup>3</sup> + 1 | *x*<sup>12647431389</sup> + *x*<sup>8008354628</sup> + *x*<sup>1277498415</sup> + 1, *x*<sup>12671067191</sup> + *x*<sup>6590666154</sup> + *x*<sup>1306461382</sup> + 1
