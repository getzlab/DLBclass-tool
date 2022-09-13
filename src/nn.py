import torch
import torch.nn as nn
import torch.nn.functional as F

class Net(nn.Module):

    def __init__(self, hlSize, inputSize, outputSize):
        super(Net, self).__init__()
        self.inputLayer = nn.Linear(inputSize, hlSize, bias=False)
        self.hl1 = nn.Linear(hlSize, hlSize)
        self.bias = nn.Parameter(torch.ones(1))
        self.outputLayer = nn.Linear(hlSize, outputSize)

    def forward(self, x):
        m = nn.Tanh()
        x = self.inputLayer(x)
        x = m(x)
        x = self.hl1(x) + self.bias
        x = m(x)
        x = self.outputLayer(x)
        x = F.softmax(x)
        return x

