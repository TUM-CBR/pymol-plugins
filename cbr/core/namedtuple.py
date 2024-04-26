from typing import Any, Callable, cast, Dict, List, NamedTuple, Optional, Tuple, Type, TypeVar

TModel = TypeVar('TModel', bound=NamedTuple)

PRIMITIVES: Dict[Type[Any], Callable[[Any], Any]] = {
    int: int,
    str: str,
    float: float
}

def get_targs(type: Type[Any]) -> Optional[Tuple[Type[Any], ...]]:
    return getattr(type, "__args__", None)

def get_optional_arg(type: Type[Any]) -> Optional[Type[Any]]:

    targs: Optional[Tuple[Type[Any], ...]] = get_targs(type)

    if targs is None or len(targs) < 1:
        return None
    
    targ = targs[0]

    if type == Optional[targ]:
        return targ
    else:
        return None
    
def get_list_arg(type: Type[Any]) -> Optional[Type[Any]]:

    targs = get_targs(type)

    if targs is None or len(targs) < 1:
        return None
    
    targ = targs[0]

    if type == List[targ]:
        return targ
    else:
        return None
    
def get_dict_args(test_type: Type[Any]) -> Optional[Tuple[Type[Any], Type[Any]]]:

    targs = get_targs(test_type)

    if targs is None or len(targs) != 2:
        return None
    
    targ0, targ1 = targs

    if test_type == Dict[targ0, targ1]:
        return targs
    else:
        return None
    
def is_namedtuple(test_type: Type[Any]) -> bool:
    return isinstance(test_type, type) \
        and issubclass(test_type, tuple) \
        and hasattr(cast(Any, test_type), "__annotations__")

def handle_value(type: Type[Any], value: Any) -> Any:

    o_type = get_optional_arg(type)

    if o_type is not None and value is not None:
        return handle_value(o_type, value)
    elif o_type is not None and value is None:
        return None

    primitive = PRIMITIVES.get(type)
    list_type = get_list_arg(type)
    dict_types = get_dict_args(type)

    if primitive is not None:
        return primitive(value)
    elif list_type is not None and isinstance(value, list):
        return [handle_value(list_type, cast(Any, v)) for v in value]
    elif dict_types is not None and isinstance(value, dict):
        t_key, t_value = dict_types
        return {
            handle_value(t_key, k): handle_value(t_value, v)
            for k,v in value.items()
        }
    elif is_namedtuple(type) and isinstance(value, dict):
        return parse(cast(Any, type), cast(Any, value))
    else:
        raise ValueError(f"Cannot parse type {type} with value {value}")

def parse(model: Type[TModel], record: Dict[Any, Any]) -> TModel:
    kwargs = {}

    for field, field_type in model.__annotations__.items():

        value = record.get(field)
        kwargs[field] = handle_value(field_type, value)

    ctor = cast(Any, model)

    return ctor(**kwargs)
