function [mu, Sigma]=kalman(mu,Sigma,obs,Q,R,A,C)
%mu��ֵ��SigmaЭ�������obs�۲�ֵ��Q����������R�۲�����
%A�����˶�ģ��ȷ���ģ�C���ɹ۲�ģ��ȷ����
[mu_bar, Sigma_bar] = kalmanPredict(mu, Sigma, A, Q);
[mu, Sigma] = kalmanUpdate(mu_bar, Sigma_bar, obs, C, R);
