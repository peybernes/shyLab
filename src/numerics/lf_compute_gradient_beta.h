      RealType beta_ooo  = beta[cell_ooo];  
      RealType beta_m1m1 = beta[cell_m1m1];
      RealType beta_p1p1 = beta[cell_p1p1];
      RealType beta_m1p1 = beta[cell_m1p1];
      RealType beta_p1m1 = beta[cell_p1m1];
      RealType beta_om1  = beta[cell_om1];
      RealType beta_op1  = beta[cell_op1];
      RealType beta_m1o  = beta[cell_m1o];
      RealType beta_p1o  = beta[cell_p1o];



      gradx = twelth * (beta_p1p1 - beta_m1p1) +
	      third  * (beta_p1o  - beta_m1o ) +
	      twelth * (beta_p1m1 - beta_m1m1);

      grady = twelth * (beta_p1p1 - beta_p1m1) +
	      third  * (beta_op1  - beta_om1 ) +
	      twelth * (beta_m1p1 - beta_m1m1);


      grad_tmp = 0.5 * grady + 0.5*gradx;
      if(grad_tmp > 0)
	{
	  zpp     = std::max(beta_ooo,beta_p1o);
	  zpp_tmp = std::max(beta_op1,beta_p1p1);
	  zpp = std::max(zpp,zpp_tmp);
	  //zsup=1.;
	  zpp = 0.5 * (zpp - beta_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpp     = std::min(beta_ooo,beta_p1o);
	  zpp_tmp = std::min(beta_op1,beta_p1p1);
	  zpp=std::min(zpp,zpp_tmp);
      //zsup=0;
	  zpp = 0.5 * (zpp - beta_ooo) / (grad_tmp);
	}
      grad_tmp = -0.5 * grady + 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zpm     = std::max(beta_ooo,beta_p1o);
	  zpm_tmp =std::max(beta_om1,beta_p1m1);
	  zpm=std::max(zpm,zpm_tmp);
	  //zsup=1.;
	  zpm = 0.5 * (zpm - beta_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zpm     = std::min(beta_ooo,beta_p1o);
	  zpm_tmp = std::min(beta_om1,beta_p1m1);
	  zpm = std::min(zpm,zpm_tmp);
      //zsup=0;
	  zpm = 0.5 * (zpm - beta_ooo) / (grad_tmp);
	}

      grad_tmp = 0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmp     = std::max(beta_ooo,beta_m1o);
	  zmp_tmp = std::max(beta_op1,beta_m1p1);
	  zmp = std::max(zmp,zmp_tmp);
	  //zsup=1.;
	  zmp = 0.5 * (zmp - beta_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmp     = std::min(beta_ooo,beta_m1o);
	  zmp_tmp = std::min(beta_op1,beta_m1p1);
	  zmp = std::min(zmp,zmp_tmp);
      //zsup=0;
	  zmp = 0.5 * (zmp - beta_ooo) / (grad_tmp);
	}

      grad_tmp = -0.5 * grady - 0.5 * gradx;
      if(grad_tmp > 0)
	{
	  zmm     = std::max(beta_ooo,beta_m1o);
	  zmm_tmp = std::max(beta_om1,beta_m1m1);
	  zmm = std::max(zmm,zmm_tmp);
	  //zsup=1.;
	  zmm = 0.5 * (zmm - beta_ooo) / (grad_tmp);
	}
      else if(grad_tmp < 0)
	{
	  zmm     = std::min(beta_ooo,beta_m1o);
	  zmm_tmp = std::min(beta_om1,beta_m1m1);
	  zmm = std::min(zmm,zmm_tmp);
      //zsup=0;
	  zmm = 0.5 * (zmm - beta_ooo) / (grad_tmp);
	}

      tmp = std::min(zpp,zmm);
      phi = std::min(zmp,zpm);
      phi = std::min(phi,tmp);
      phi = std::min(phi,2.);
