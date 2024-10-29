import torch
from torch_geometric.typing import OptTensor, PairOptTensor, PairTensor, Adj
from torch_sparse import SparseTensor, set_diag
from torch_geometric.utils import remove_self_loops, add_self_loops
from torch_geometric.nn.conv import MessagePassing
from typing import Optional, Callable, Union
from torch import Tensor
# from ..inits import
from torch_geometric.nn.inits import  reset

def get_angle(v1: Tensor, v2: Tensor) -> Tensor:
    return torch.atan2(
        torch.cross(v1, v2, dim=1).norm(p=2, dim=1), (v1 * v2).sum(dim=1))


def point_pair_features(pos_i: Tensor, pos_j: Tensor, normal_i: Tensor,
                        normal_j: Tensor) -> Tensor:
    pseudo = pos_j - pos_i
    return torch.stack([
        pseudo.norm(p=2, dim=1),
        get_angle(normal_i, pseudo),
        get_angle(normal_j, pseudo),
        get_angle(normal_i, normal_j)
    ], dim=1)


class PPFConv(MessagePassing):
    def __init__(self, local_nn: Optional[Callable] = None,
                 global_nn: Optional[Callable] = None,
                 add_self_loops: bool = True, **kwargs):
        kwargs.setdefault('aggr', 'max')
        super(PPFConv, self).__init__(**kwargs)

        self.local_nn = local_nn
        self.global_nn = global_nn
        self.add_self_loops = add_self_loops

        self.reset_parameters()

    def reset_parameters(self):
        reset(self.local_nn)
        reset(self.global_nn)

    def normalized_point_pair_features(self, pos_i: Tensor, pos_j: Tensor, normal_i: Tensor,
                            normal_j: Tensor, radius : float) -> Tensor:
        pseudo = pos_j - pos_i
        angle1 = get_angle(normal_i, pseudo)
        angle2 = get_angle(normal_j, pseudo)
        angle3 = get_angle(normal_i, normal_j)
        return torch.stack([
            # pseudo.norm(p=2, dim=1) / radius,
            pseudo.norm(p=2, dim=1) ,
            torch.sin(angle1),
            # torch.cos(angle1),
            torch.sin(angle2),
            # torch.cos(angle2),
            torch.sin(angle3),
            # torch.cos(angle3)
        ], dim=1)

    def forward(self, x: Union[OptTensor, PairOptTensor],
                pos: Union[Tensor, PairTensor],
                normal: Union[Tensor, PairTensor],
                edge_index: Adj,
                radius: float) -> Tensor:  # yapf: disable
        """"""
        if not isinstance(x, tuple):
            x: PairOptTensor = (x, None)

        if isinstance(pos, Tensor):
            pos: PairTensor = (pos, pos)

        if isinstance(normal, Tensor):
            normal: PairTensor = (normal, normal)

        if self.add_self_loops:
            if isinstance(edge_index, Tensor):
                edge_index, _ = remove_self_loops(edge_index)
                edge_index, _ = add_self_loops(edge_index,
                                               num_nodes=pos[1].size(0))
            elif isinstance(edge_index, SparseTensor):
                edge_index = set_diag(edge_index)

        # propagate_type: (x: PairOptTensor, pos: PairTensor, normal: PairTensor)  # noqa
        out = self.propagate(edge_index, x=x, pos=pos, normal=normal,
                             size=None, radius= radius)

        if self.global_nn is not None:
            out = self.global_nn(out)

        return out

    def message(self, x_j: OptTensor, pos_i: Tensor, pos_j: Tensor,
                normal_i: Tensor, normal_j: Tensor, radius_i : float) -> Tensor:
        msg = self.normalized_point_pair_features(pos_i, pos_j, normal_i, normal_j, radius_i)
        if x_j is not None:
            msg = torch.cat([x_j, msg], dim=1)
        if self.local_nn is not None:
            msg = self.local_nn(msg)
        return msg

    def __repr__(self):
        return '{}(local_nn={}, global_nn={})'.format(self.__class__.__name__,
                                                      self.local_nn,
                                                      self.global_nn)

