import torch.nn as nn
from collections import OrderedDict

# define NN architecture
class twolayermlp(nn.Module):
    def __init__(self,**kwargs):
        super(twolayermlp, self).__init__()
        in_channels = kwargs.get('num_channels', 85)
        classes = kwargs.get('num_classes', 2)

        self.model = nn.Sequential(OrderedDict([
            ('fc1', nn.Linear(in_channels, 512)),
            ('relu1', nn.ReLU()),
            ('dropout1',nn.Dropout(0.1)),
            ('fc2', nn.Linear(512, 10)),
            ('relu2', nn.ReLU()),
            ('dropout2', nn.Dropout(0.1)),
            ('output', nn.Linear(10, classes)),
            ('softmax', nn.Softmax(dim=1))]))


    def forward(self, x):
        output=self.model(x)
        return output


def mlp(**kwargs):
    model = twolayermlp(**kwargs)
    return model

