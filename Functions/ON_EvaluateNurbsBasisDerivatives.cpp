bool ON_EvaluateNurbsBasisDerivatives(
  int order,
  const double* knot, 
  int der_count,
  double* N 
)
{
	double dN, c;
	const double *k0, *k1;
	double *a0, *a1, *ptr, **dk;
	int i, j, k, jmax;

	const int d = order-1;
	const int Nstride = -der_count*order;

	/* workspaces for knot differences and coefficients 
	 *
	 * a0[] and a1[] have order doubles
	 *
	 * dk[0] = array of d knot differences
	 * dk[1] = array of (d-1) knot differences
	 *
	 * dk[der_count-1] = 1.0/(knot[d] - knot[d-1])
	 * dk[der_count] = dummy pointer to make loop efficient
	 */

	//dk = (double**)alloca( (der_count+1) << 3 ); /* << 3 in case pointers are 8 bytes long */
	//a0 = (double*)alloca( (order*(2 + ((d+1)>>1))) << 3 ); /* d for a0, d for a1, d*order/2 for dk[]'s and slop to avoid /2 */
	//a1 = a0 + order;

  double stack_buffer[80];
  void* heap_buffer = 0;
  const size_t dbl_count = (order*(2 + ((d+1)>>1)));
  const size_t sz = ( dbl_count*sizeof(*a0) + (der_count+1)*sizeof(*dk) );

  a0 = (sz <= sizeof(stack_buffer)) ? stack_buffer : (double*)(heap_buffer = onmalloc(sz));
	a1 = a0 + order;
  dk = (double**)(a0 + dbl_count);

	/* initialize reciprocal of knot differences */
	dk[0] = a1 + order;
	for (k = 0; k < der_count; k++) {
		j = d-k;
		k0 = knot++;
		k1 = k0 + j;
		for (i = 0; i < j; i++) 
			dk[k][i] = 1.0/(*k1++ - *k0++);
		dk[k+1] = dk[k] + j;
	}
	dk--;
	/* dk[1] = 1/{t[d]-t[0], t[d+1]-t[1], ..., t[2d-2] - t[d-2], t[2d-1] - t[d-1]}
	 *       = diffs needed for 1rst derivative
	 * dk[2] = 1/{t[d]-t[1], t[d+1]-t[2], ..., t[2d-2] - t[d-1]}
	 *       = diffs needed for 2nd derivative
	 * ...
	 * dk[d] = 1/{t[d] - t[d-1]}
	 *       = diff needed for d-th derivative
	 *
	 * d[k][n] = 1.0/( t[d+n] - t[k-1+n] )
	 */
	
	N += order;
	/* set N[0] ,..., N[d] = 1rst derivatives, 
	 * N[order], ..., N[order+d] = 2nd, etc.
	 */
	for ( i=0; i<order; i++) {
		a0[0] = 1.0;
		for (k = 1; k <= der_count; k++) {
			/* compute k-th derivative of N_i^d up to d!/(d-k)! scaling factor */
			dN = 0.0;
			j = k-i; 
			if (j <= 0) {
				dN = (a1[0] = a0[0]*dk[k][i-k])*N[i];
				j = 1;
			}
			jmax = d-i; 
			if (jmax < k) {
				while (j <= jmax) {
					dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
					j++;
				}
			}
			else {
				/* sum j all the way to j = k */
				while (j < k) {
					dN += (a1[j] = (a0[j] - a0[j-1])*dk[k][i+j-k])*N[i+j];
					j++;
				}
				dN += (a1[k] = -a0[k-1]*dk[k][i])*N[i+k];
			}

			/* d!/(d-k)!*dN = value of k-th derivative */
			N[i] = dN;
			N += order;
			/* a1[] s for next derivative = linear combination
			 * of a[]s used to compute this derivative.
			 */
			ptr = a0; a0 = a1; a1 = ptr;
		}
		N += Nstride;
	}

	/* apply d!/(d-k)! scaling factor */
	dN = c = (double)d;
	k = der_count;
	while (k--) {
		i = order;
		while (i--)
			*N++ *= c;
		dN -= 1.0;
		c *= dN;
	}

  if ( 0 != heap_buffer )
    onfree(heap_buffer);

  return true;
}
