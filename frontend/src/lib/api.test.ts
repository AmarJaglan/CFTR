import {describe,expect,it,vi} from 'vitest';
import {api} from './api';

describe('typed API client',()=>{it('encodes variant identifiers',async()=>{const fetchMock=vi.spyOn(globalThis,'fetch').mockResolvedValue(new Response(JSON.stringify({canonical_protein:'F508del'}),{status:200}));await api.getVariant('F 508/del');expect(fetchMock.mock.calls[0][0]).toContain('/variants/F%20508%2Fdel');fetchMock.mockRestore()});it('surfaces backend errors',async()=>{vi.spyOn(globalThis,'fetch').mockResolvedValue(new Response(JSON.stringify({detail:'not found'}),{status:404}));await expect(api.getVariant('nope')).rejects.toThrow('not found');vi.restoreAllMocks()})});
