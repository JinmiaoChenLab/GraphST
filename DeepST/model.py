import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
from torch_geometric.nn import GCNConv

class Discriminator(nn.Module):
    def __init__(self, n_h):
        super(Discriminator, self).__init__()
        self.f_k = nn.Bilinear(n_h, n_h, 1)

        for m in self.modules():
            self.weights_init(m)

    def weights_init(self, m):
        if isinstance(m, nn.Bilinear):
            torch.nn.init.xavier_uniform_(m.weight.data)
            if m.bias is not None:
                m.bias.data.fill_(0.0)

    def forward(self, c, h_pl, h_mi, s_bias1=None, s_bias2=None):
        c_x = c.expand_as(h_pl)  

        sc_1 = self.f_k(h_pl, c_x)
        sc_2 = self.f_k(h_mi, c_x)

        if s_bias1 is not None:
            sc_1 += s_bias1
        if s_bias2 is not None:
            sc_2 += s_bias2

        logits = torch.cat((sc_1, sc_2), 1)

        return logits
    
class AvgReadout(nn.Module):
    def __init__(self):
        super(AvgReadout, self).__init__()

    def forward(self, seq, msk=None):
        if msk is None:
            return torch.mean(seq, 0)
        else:
            msk = torch.unsqueeze(msk, -1)
            return torch.sum(seq * msk, 0) / torch.sum(msk) 
    
class Encoder(Module):
    """
    Simple GCN layer, similar to https://arxiv.org/abs/1609.02907
    """
    def __init__(self, in_features, out_features, dropout=0.0, act=F.relu):
        super(Encoder, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.dropout = dropout
        self.act = act
        
        self.weight1 = Parameter(torch.FloatTensor(self.in_features, self.out_features))
        self.weight2 = Parameter(torch.FloatTensor(self.out_features, self.in_features))
        self.reset_parameters()
        
        self.disc = Discriminator(out_features)

        self.sigm = nn.Sigmoid()
        self.read = AvgReadout()
        
    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.weight1)
        torch.nn.init.xavier_uniform_(self.weight2)

    def forward(self, input, input_a, adj):
        input = F.dropout(input, self.dropout, self.training)
        support1 = torch.mm(input, self.weight1)
        output1 = torch.mm(adj, support1)
        
        hiden_feat = output1
        
        support2 = torch.mm(output1, self.weight2)
        emb = torch.mm(adj, support2)
        
        output1 = self.act(output1)
        
        #adversial learning
        input_a = F.dropout(input_a, self.dropout, self.training)
        support_a1 = torch.mm(input_a, self.weight1)
        output_a1 = torch.mm(adj, support_a1)
        output_a1 = self.act(output_a1)
        
        h = self.read(output1) 
        h = self.sigm(h)  

        h_a = self.read(output_a1)
        h_a = self.sigm(h_a)  

        # adversarial learning
        ret = self.disc(h, output1, output_a1)  
        ret_a = self.disc(h_a, output_a1, output1) 
        
        return hiden_feat, emb, ret, ret_a 

class Encoder_sc(torch.nn.Module):
    def __init__(self, dim_input, dim_output, act=F.elu):
        super(Encoder_sc, self).__init__()
        self.dim_input = dim_input
        self.dim_output = dim_output
        self.act = act
        
        self.weight1_en = Parameter(torch.FloatTensor(self.dim_input, self.dim_output))
        self.weight1_de = Parameter(torch.FloatTensor(self.dim_output, self.dim_input))
  
        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.weight1_en)
        torch.nn.init.xavier_uniform_(self.weight1_de)
        
    def forward(self, x):
        x = torch.mm(x, self.weight1_en)
        x = torch.mm(x, self.weight1_de)
        
        return x
    
class Encoder_map(torch.nn.Module):
    def __init__(self, args, n_cell, n_spot, n_domain=100):
        super(Encoder_map, self).__init__()
        self.n_cell = n_cell
        self.n_spot = n_spot
        self.n_domain = n_domain
        
        if args.domain:
           self.M = Parameter(torch.FloatTensor(self.n_cell, self.n_domain)) 
        else:    
           self.M = Parameter(torch.FloatTensor(self.n_cell, self.n_spot))
        
        self.reset_parameters()

    def reset_parameters(self):
        torch.nn.init.xavier_uniform_(self.M)
        
    def forward(self):
        x = self.M
        
        return x 
   
    
   
    

    
    
        
       
