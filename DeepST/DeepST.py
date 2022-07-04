import torch
from preprocess import preprocess_adj, preprocess, construct_interaction, add_contrastive_label, get_feature
import time
import random
import numpy as np
from model import Encoder, Encoder_map, Encoder_sc
from progress.bar import Bar
from torch import nn
from preprocess import permutation, fix_seed
import torch.nn.functional as F
from scipy.sparse.csc import csc_matrix
from scipy.sparse.csr import csr_matrix
import pandas as pd

class DeepST():
    def __init__(self, adata, random_seed=50, add_regularization=True):
       self.adata = adata.copy()
       self.random_seed = random_seed
       self.add_regularization = True
       
       fix_seed(self.random_seed)
       preprocess(self.adata)
       construct_interaction(self.adata)
       add_contrastive_label(self.adata)
       
       self.adata_output = self.adata.copy()
    
    def train_DeepST(self):
       if self.add_regularization:
          adata = self.adata_output.copy()
          #preprocess(adata)
          get_feature(adata)
          model = Train(adata)
          emb = model.train()
          self.adata_output.obsm['emb'] = emb
          
          fix_seed(self.random_seed)
          adata = self.adata_output.copy()
          #preprocess(adata)
          get_feature(adata)
          model = Train(adata, add_regularization=True)
          emb_regularization = model.train()
          self.adata_output.obsm['emb_reg'] = emb_regularization
          
       else:
          model = Train(self.adata.copy())
          emb = model.train()
          self.adata_output.obsm['emb'] = emb
          
       return self.adata_output   
    
class Train():
    def __init__(self, 
            adata,
            adata_sc = None,
            device='cuda:0',
            learning_rate=0.001,
            learning_rate_sc = 0.01,
            weight_decay=0.00,
            epochs=600, 
            dim_input=3000,
            dim_output=64,
            random_seed = 50,
            alpha = 10,
            beta = 1,
            theta = 0.1,
            lamda1 = 10,
            lamda2 = 1,
            add_regularization = False,
            deconvolution = False
            ):
        self.adata = adata.copy()
        self.device = device
        self.learning_rate=learning_rate
        self.learning_rate_sc = learning_rate_sc
        self.weight_decay=weight_decay
        self.epochs=epochs
        self.alpha = alpha
        self.beta = beta
        self.theta = theta
        self.lamda1 = lamda1
        self.lamda2 = lamda2
        self.add_regularization = add_regularization
        self.deconvolution = deconvolution
        
        self.features = torch.FloatTensor(adata.obsm['feat'].copy()).to(self.device)
        self.features_a = torch.FloatTensor(adata.obsm['feat_a'].copy()).to(self.device)
        self.label_CSL = torch.FloatTensor(adata.obsm['label_CSL']).to(self.device)
        self.adj = adata.obsm['adj']
        self.graph_neigh = torch.FloatTensor(adata.obsm['graph_neigh'].copy() + np.eye(self.adj.shape[0])).to(self.device)
    
        self.dim_input = self.features.shape[1]
        self.dim_output = dim_output
        
        # standard version
        self.adj = preprocess_adj(self.adj)
        self.adj = torch.FloatTensor(self.adj).to(self.device)
        
        if self.deconvolution:
           self.adata_sc = adata_sc.copy() 
            
           if isinstance(self.adata.X, csc_matrix) or isinstance(self.adata.X, csr_matrix):
              self.feat_sp = adata.X.toarray()[:, ]
           else:
              self.feat_sp = adata.X[:, ]
           if isinstance(self.adata_sc.X, csc_matrix) or isinstance(self.adata_sc.X, csr_matrix):
              self.feat_sc = self.adata_sc.X.toarray()[:, ]
           else:
              self.feat_sc = self.adata_sc.X[:, ]
            
           # fill nan as 0
           self.feat_sc = pd.DataFrame(self.feat_sc).fillna(0).values
           self.feat_sp = pd.DataFrame(self.feat_sp).fillna(0).values
          
           self.feat_sc = torch.FloatTensor(self.feat_sc).to(self.device)
           self.feat_sp = torch.FloatTensor(self.feat_sp).to(self.device)
        
           if self.adata_sc is not None:
              self.dim_input = self.feat_sc.shape[1] 

           self.n_cell = adata_sc.n_obs
           self.n_spot = adata.n_obs
            
    def train(self):
        self.model = Encoder(self.dim_input, self.dim_output, self.graph_neigh).to(self.device)
        self.loss_CSL = nn.BCEWithLogitsLoss()
    
        self.optimizer = torch.optim.Adam(self.model.parameters(), self.learning_rate, 
                                          weight_decay=self.weight_decay)
        
        t = time.time()
        if not self.add_regularization:
           print('Begin to train ST data.')
        bar = Bar(max=self.epochs)
        bar.check_tty = False
        self.model.train()
       
        for epoch in range(self.epochs):
            self.model.train()
              
            self.features_a = permutation(self.features)
            self.hiden_feat, self.emb, ret, ret_a = self.model(self.features, self.features_a, self.adj)
            
            self.loss_sl_1 = self.loss_CSL(ret, self.label_CSL)
            self.loss_sl_2 = self.loss_CSL(ret_a, self.label_CSL)
            self.loss_feat = F.mse_loss(self.features, self.emb)
            
            if self.add_regularization:
               self.loss_norm = 0
               for name, parameters in self.model.named_parameters():
                   if name in ['weight1', 'weight2']:
                      self.loss_norm = self.loss_norm + torch.norm(parameters, p=2) 
               loss =  self.alpha*self.loss_feat + self.beta*(self.loss_sl_1 + self.loss_sl_2) + self.theta*self.loss_norm 
            else: 
               loss =  self.alpha*self.loss_feat + self.beta*(self.loss_sl_1 + self.loss_sl_2)
            
            self.optimizer.zero_grad()
            loss.backward() 
            self.optimizer.step()
            
            if not self.add_regularization:
               bar_str = '{}/{}|Loss:{:.10f}|Time:{:.4f}'
               bar.suffix = bar_str.format(epoch + 1, self.epochs, 
                                        loss.detach().cpu().numpy(), time.time() - t)
               bar.next()
                      
        bar.finish()
        
        if not self.add_regularization:
           print("Training finished for spot representation learning!")
        
        with torch.no_grad():
             self.model.eval()
             if self.deconvolution:
                self.emb_rec = self.model(self.features, self.features_a, self.adj)[1]
             else:   
                self.emb_rec = self.model(self.features, self.features_a, self.adj)[1].detach().cpu().numpy()
             
             return self.emb_rec
         
    def train_sc(self):
        self.model_sc = Encoder_sc(self.dim_input, self.dim_output).to(self.device)
        self.optimizer_sc = torch.optim.Adam(self.model_sc.parameters(), lr=self.learning_rate_sc)  
        
        print('Begin to train scRNA data.')
        bar = Bar(max=self.epochs)
        bar.check_tty = False
        t = time.time()
        for epoch in range(self.epochs):
            self.model_sc.train()
            
            emb = self.model_sc(self.feat_sc)
            loss = F.mse_loss(emb, self.feat_sc)
            
            self.optimizer_sc.zero_grad()
            loss.backward()
            self.optimizer_sc.step()
            
            bar_str = '{}/{}|Loss:{:.10f}|Time:{:.4f}'
            bar.suffix = bar_str.format(epoch + 1, self.epochs, 
                                        loss.detach().cpu().numpy(), time.time() - t)
            bar.next()
        bar.finish()
        print("Optimization finished for cell representation learning!")
        
        with torch.no_grad():
            self.model_sc.eval()
            emb_sc = self.model_sc(self.feat_sc)
         
            return emb_sc
        
    def train_map(self):
        emb_sp = self.train()
        emb_sc = self.train_sc()
        
        self.adata.obsm['emb_sp'] = emb_sp.detach().cpu().numpy()
        self.adata_sc.obsm['emb_sc'] = emb_sc.detach().cpu().numpy()
        
        # Normalize features for consistence between ST and scRNA-seq
        emb_sp = F.normalize(emb_sp, p=2, eps=1e-12, dim=1)
        emb_sc = F.normalize(emb_sc, p=2, eps=1e-12, dim=1)
        
        self.model_map = Encoder_map(self.n_cell, self.n_spot).to(self.device)  
          
        self.optimizer_map = torch.optim.Adam(self.model_map.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
        
        print('Begin to learn mapping matrix.')
        bar = Bar(max=self.epochs)
        bar.check_tty = False
        t = time.time()
        for epoch in range(self.epochs):
            self.model_map.train()
            self.map_matrix = self.model_map()

            loss_recon, loss_NCE = self.loss(emb_sp, emb_sc)
             
            loss = self.lamda1*loss_recon + self.lamda2*loss_NCE 

            self.optimizer_map.zero_grad()
            loss.backward()
            self.optimizer_map.step()
            
            bar_str = '{}/{}|Loss:{:.10f}|Time:{:.4f}'
            bar.suffix = bar_str.format(epoch + 1, self.epochs, 
                                        loss.detach().cpu().numpy(), time.time() - t)
            bar.next()
        bar.finish()
        print("Mapping matrix learning finished!")
        
        # take final softmax w/o computing gradients
        with torch.no_grad():
            self.model_map.eval()
            emb_sp = emb_sp.cpu().numpy()
            emb_sc = emb_sc.cpu().numpy()
            map_matrix = F.softmax(self.map_matrix, dim=0).cpu().numpy() # dim=0: normalization by spot
            
            self.adata.obsm['emb_sp'] = emb_sp
            self.adata_sc.obsm['emb_sc'] = emb_sc
            self.adata.obsm['map_matrix'] = map_matrix.T # spot x cell

            return self.adata, self.adata_sc
    
    def loss(self, emb_sp, emb_sc):
        # cell-to-spot
        map_probs = F.softmax(self.map_matrix, dim=0)   # dim=0: normalization by spot
        self.pred_sp = torch.matmul(map_probs.t(), emb_sc)
           
        loss_recon = F.mse_loss(self.pred_sp, emb_sp, reduction='mean')
        loss_NCE = self.Noise_Cross_Entropy(self.pred_sp, emb_sp)
           
        return loss_recon, loss_NCE
        
    def Noise_Cross_Entropy(self, pred_sp, emb_sp):
        '''
        Considering spatial neighbors as positive pairs for each spot
        '''
        mat = self.cosine_similarity(pred_sp, emb_sp) 
        k = torch.exp(mat).sum(axis=1) - torch.exp(torch.diag(mat, 0))
        
        # positive pairs
        p = torch.exp(mat)
        p = torch.mul(p, self.graph_neigh).sum(axis=1)
        
        ave = torch.div(p, k)
        loss = - torch.log(ave).mean()
        
        return loss
    
    def cosine_similarity(self, pred_sp, emb_sp):  #pres_sp: spot x gene; emb_sp: spot x gene    
        M = torch.matmul(pred_sp, emb_sp.T)
        Norm_c = torch.norm(pred_sp, p=2, dim=1)
        Norm_s = torch.norm(emb_sp, p=2, dim=1)
        Norm = torch.matmul(Norm_c.reshape((pred_sp.shape[0], 1)), Norm_s.reshape((emb_sp.shape[0], 1)).T) + -5e-12
        M = torch.div(M, Norm)
        
        if torch.any(torch.isnan(M)):
           M = torch.where(torch.isnan(M), torch.full_like(M, 0.4868), M)

        return M        