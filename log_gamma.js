var gammaCof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
                -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];


function log_gamma(xx)
{
    var x = xx - 1.0;
    var tmp = x + 5.5; tmp -= (x + 0.5)*Math.log(tmp);
    var ser=1.000000000190015;
    for (j=0;j<=5;j++){ x++; ser += gammaCof[j]/x; }
    return -tmp+Math.log(2.5066282746310005*ser);
}


function dirichletScore(params,vals){
	params.pop()
	vals.pop()
	var alpha = params;
	var theta = vals;

	var asum = 0;
	for (var i = 0; i < alpha.length; i++) {
	  asum += alpha[i];
	}
	var logp = log_gamma(asum);

	for (var i = 0; i < alpha.length; i++){

	   logp += (alpha[i]-1)*Math.log(theta[i]);
	   logp -= log_gamma(alpha[i]);
	 }
	 return logp;
}