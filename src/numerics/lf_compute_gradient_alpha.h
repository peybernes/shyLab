      sum_alphax = 0;
      sum_alphay = 0;

      for (int k = 0;k<nb_mat;k++) {

	RealType alpha_ooo;  
	RealType alpha_m1m1;
	RealType alpha_p1p1;
	RealType alpha_m1p1;
	RealType alpha_p1m1;
	RealType alpha_om1;
	RealType alpha_op1;
	RealType alpha_m1o;
	RealType alpha_p1o;
	
	if (k==0) {
	  alpha_ooo  = in_c_1[cell_ooo] ;
	  alpha_m1m1 = in_c_1[cell_m1m1];
	  alpha_p1p1 = in_c_1[cell_p1p1];
	  alpha_m1p1 = in_c_1[cell_m1p1];
	  alpha_p1m1 = in_c_1[cell_p1m1];
	  alpha_om1  = in_c_1[cell_om1] ;
	  alpha_op1  = in_c_1[cell_op1] ;
	  alpha_m1o  = in_c_1[cell_m1o] ;
	  alpha_p1o  = in_c_1[cell_p1o] ;
	}
	else {
	  alpha_ooo  = in_c_2[cell_ooo] ;
	  alpha_m1m1 = in_c_2[cell_m1m1];
	  alpha_p1p1 = in_c_2[cell_p1p1];
	  alpha_m1p1 = in_c_2[cell_m1p1];
	  alpha_p1m1 = in_c_2[cell_p1m1];
	  alpha_om1  = in_c_2[cell_om1] ;
	  alpha_op1  = in_c_2[cell_op1] ;
	  alpha_m1o  = in_c_2[cell_m1o] ;
	  alpha_p1o  = in_c_2[cell_p1o] ;
	}
	
	zpp=2.;
	zpm=2.;
	zmp=2.;
	zmm=2.;
	
	gradx = twelth * (alpha_p1p1 - alpha_m1p1 ) +
	  third  * (alpha_p1o  - alpha_m1o  ) +
	  twelth * (alpha_p1m1 - alpha_m1m1 );
	
	grady = twelth * (alpha_p1p1 - alpha_p1m1 ) +
	  third  * (alpha_op1  - alpha_om1  ) +
	  twelth * (alpha_m1p1 - alpha_m1m1 );
	
	grad_tmp = 0.5 * (grady + gradx);
	
	if(grad_tmp > 0)
	  {
	    zpp     = std::max(alpha_ooo,alpha_p1o);
	    zpp_tmp = std::max(alpha_op1,alpha_p1p1);
	    zpp     = std::max(zpp,zpp_tmp);
	    zpp     = (zpp - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zpp     = std::min(alpha_ooo,alpha_p1o);
	    zpp_tmp = std::min(alpha_op1,alpha_p1p1);
	    zpp     = std::min(zpp,zpp_tmp);
	    zpp     = (zpp - alpha_ooo) / (grad_tmp);
	  }
	
	grad_tmp = 0.5 * (-grady + gradx);
	
	if(grad_tmp > 0)
	  {
	    zpm     = std::max(alpha_ooo,alpha_p1o);
	    zpm_tmp = std::max(alpha_om1,alpha_p1m1);
	    zpm     = std::max(zpm,zpm_tmp);
	    zpm     = (zpm - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zpm     = std::min(alpha_ooo,alpha_p1o);
	    zpm_tmp = std::min(alpha_om1,alpha_p1m1);
	    zpm     = std::min(zpm,zpm_tmp);
	    zpm     = (zpm - alpha_ooo) / (grad_tmp);
	  }
	
	grad_tmp = 0.5 * (grady - gradx);
	
	if(grad_tmp > 0)
	  {
	    zmp     = std::max(alpha_ooo,alpha_m1o);
	    zmp_tmp = std::max(alpha_op1,alpha_m1p1);
	    zmp     = std::max(zmp,zmp_tmp);
	    zmp     = (zmp - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zmp     = std::min(alpha_ooo,alpha_m1o);
	    zmp_tmp = std::min(alpha_op1,alpha_m1p1);
	    zmp     = std::min(zmp,zmp_tmp);
	    zmp     = (zmp - alpha_ooo) / (grad_tmp);
	  }
	
	grad_tmp = -0.5 * (grady + gradx);
	
	if(grad_tmp > 0)
	  {
	    zmm     = std::max(alpha_ooo,alpha_m1o);
	    zmm_tmp = std::max(alpha_om1,alpha_m1m1);
	    zmm     = std::max(zmm,zmm_tmp);
	    zmm     = (zmm - alpha_ooo) / (grad_tmp);
	  }
	else if(grad_tmp < 0)
	  {
	    zmm     = std::min(alpha_ooo,alpha_m1o);
	    zmm_tmp = std::min(alpha_om1,alpha_m1m1);
	    zmm     = std::min(zmm,zmm_tmp);
	    zmm     = (zmm - alpha_ooo) / (grad_tmp);
	  }
	
	tmp = std::min(zpp,zmm);
	phi = std::min(zmp,zpm);
	phi = std::min(phi,tmp);
	phi = std::min(phi,2.);
	
	sum_alphax = sum_alphax + 0.5 * phi * gradx;
	sum_alphay = sum_alphay + 0.5 * phi * grady;
	
	alphak_grady_top  [k][cell_ooo] = phi * grady * h_y;
	alphak_gradx_right[k][cell_ooo] = phi * gradx * h_x;
      
      }

      for (int k = 0;k<nb_mat;k++) {
      
	RealType alpha_ooo;  

	if (k==0) {
	  alpha_ooo  = in_c_1[cell_ooo] ;
	}
	else {
	  alpha_ooo  = in_c_2[cell_ooo] ;
	}

	alphak_grady_bot  [k][cell_ooo] = -2 * h_y * (alpha_ooo * (sum_alphay)  - 0.5 * dx * alphak_grady_top  [k][cell_ooo]) / (1 - sum_alphay);
      
	alphak_gradx_left [k][cell_ooo] = -2 * h_x * (alpha_ooo * (sum_alphax)  - 0.5 * dx * alphak_gradx_right[k][cell_ooo]) / (1 - sum_alphax);
      
	alphak_grady_top  [k][cell_ooo] =  2 * h_y * (alpha_ooo * (-sum_alphay) + 0.5 * dx * alphak_grady_top  [k][cell_ooo]) / (1 + sum_alphay);
      
	alphak_gradx_right[k][cell_ooo] =  2 * h_x * (alpha_ooo * (-sum_alphax) + 0.5 * dx * alphak_gradx_right[k][cell_ooo]) / (1 + sum_alphax);

      }
