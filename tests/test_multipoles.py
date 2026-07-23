import pytest

from acoustic_ms.multipoles import mode_count, mode_from_index, mode_index, modes


def test_mode_count_order_and_round_trip():
    assert mode_count(3) == 16
    assert modes(1) == ((0, 0), (1, -1), (1, 0), (1, 1))
    for index in range(mode_count(5)):
        assert mode_index(*mode_from_index(index)) == index


@pytest.mark.parametrize("value", [-1, 1.0])
def test_invalid_mode_inputs_are_rejected(value):
    with pytest.raises(ValueError):
        mode_count(value)
    with pytest.raises(ValueError):
        mode_from_index(value)
