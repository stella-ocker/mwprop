#!/usr/bin/env python

# mwprop v2.0 Jan 2026

'''

tests NE2025p runs properly and gives expected output

'''

from mwprop.nemod import *
from mwprop.nemod.NE2025 import ne2025


def main():
	# J0323+3944, D_obs = 0.95, DM_obs = 26.19
	test_gl = 152.180
	test_gb = -14.338
	test_dm_from_d = 16.2718
	test_dkpc_from_dm = 1.3863
	test_sm_from_d = 1.3794e-04
	test_sm_from_dm = 3.5591e-04

	Dk1,Dv1,Du1,Dd1 = ne2025(test_gl,test_gb,0.95,-1,classic=False,dmd_only=False)
	Dk2,Dv2,Du2,Dd2 = ne2025(test_gl,test_gb,26.19,1,classic=False,dmd_only=False)

	print('Testing NE2025p integrations give expected outputs ... Percent errors listed below (values should be less than a few percent):')

	try:
		DM1 = Dv1['DM']
		SM1 = Dv1['SM']
		DM_err = 100*abs(DM1-test_dm_from_d)/test_dm_from_d
		SM_err = 100*abs(SM1-test_sm_from_d)/test_sm_from_d
		print('D->DM: ', DM_err)
		print('D->SM: ', SM_err)

	except:
		print('Error: cannot find expected output (D->DM)')

	try:
		DIST2 = Dv2['DIST']
		SM2 = Dv2['SM']
		Derr = 100*abs(DIST2-test_dkpc_from_dm)/test_dkpc_from_dm
		SMerr = 100*abs(SM2-test_sm_from_dm)/test_sm_from_dm
		print('DM->D: ', Derr)
		print('DM->SM: ', SM_err)

	except:
		print('Error: cannot find expected output (DM->D)')


if __name__ == '__main__':
	main()